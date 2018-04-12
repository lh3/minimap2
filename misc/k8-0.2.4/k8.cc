/* The MIT License

   Copyright (c) 2012-2013, 2016 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#define K8_VERSION "0.2.4-r79" // known to work with V8-3.16.14

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <zlib.h>
#include <v8.h>

/************************************
 *** Convenient macros from v8cgi ***
 ************************************/

// A V8 object can have multiple internal fields invisible to JS. The following 4 macros read and write these fields
#define SAVE_PTR(_args, _index, _ptr)  (_args).This()->SetAlignedPointerInInternalField(_index, (void*)(_ptr))
#define LOAD_PTR(_args, _index, _type) reinterpret_cast<_type>((_args).This()->GetAlignedPointerFromInternalField(_index))
#define SAVE_VALUE(_args, _index, _val) (_args).This()->SetInternalField(_index, _val)
#define LOAD_VALUE(_args, _index) (_args).This()->GetInternalField(_index)

#define JS_STR(...) v8::String::New(__VA_ARGS__)

#define JS_THROW(type, reason) v8::ThrowException(v8::Exception::type(JS_STR(reason)))
#define JS_ERROR(reason) JS_THROW(Error, reason)
#define JS_METHOD(_func, _args) v8::Handle<v8::Value> _func(const v8::Arguments &(_args))

#define ASSERT_CONSTRUCTOR(_args) if (!(_args).IsConstructCall()) { return JS_ERROR("Invalid call format. Please use the 'new' operator."); }

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

/*******************************
 *** Fundamental V8 routines ***
 *******************************/

static inline const char *k8_cstr(const v8::String::AsciiValue &str) // Convert a V8 string to C string
{
	return *str? *str : "<N/A>";
}

static void k8_exception(v8::TryCatch *try_catch) // Exception handling. Adapted from v8/shell.cc
{
	v8::HandleScope handle_scope;
	v8::String::AsciiValue exception(try_catch->Exception());
	const char* exception_string = k8_cstr(exception);
	v8::Handle<v8::Message> message = try_catch->Message();
	if (message.IsEmpty()) {
		// V8 didn't provide any extra information about this error; just print the exception.
		fprintf(stderr, "%s\n", exception_string);
	} else {
		// Print (filename):(line number): (message).
		v8::String::AsciiValue filename(message->GetScriptResourceName());
		const char* filename_string = k8_cstr(filename);
		int linenum = message->GetLineNumber();
		fprintf(stderr, "%s:%i: %s\n", filename_string, linenum, exception_string);
		// Print line of source code.
		v8::String::AsciiValue sourceline(message->GetSourceLine());
		const char *sourceline_string = k8_cstr(sourceline);
		fprintf(stderr, "%s\n", sourceline_string);
		// Print wavy underline (GetUnderline is deprecated).
		int start = message->GetStartColumn();
		for (int i = 0; i < start; i++) fputc(' ', stderr);
		int end = message->GetEndColumn();
		for (int i = start; i < end; i++) fputc('^', stderr);
		fputc('\n', stderr);
		v8::String::AsciiValue stack_trace(try_catch->StackTrace());
		if (stack_trace.length() > 0) { // TODO: is the following redundant?
			const char* stack_trace_string = k8_cstr(stack_trace);
			fputs(stack_trace_string, stderr); fputc('\n', stderr);
		}
	}
}

bool k8_execute(v8::Handle<v8::String> source, v8::Handle<v8::Value> name, bool prt_rst) // Execute JS in a string. Adapted from v8/shell.cc
{
	v8::HandleScope handle_scope;
	v8::TryCatch try_catch;
	if (source == v8::Handle<v8::String>()) return false;
	v8::Handle<v8::Script> script = v8::Script::Compile(source, name);
	if (script.IsEmpty()) {
		k8_exception(&try_catch);
		return false;
	} else {
		v8::Handle<v8::Value> result = script->Run();
		if (result.IsEmpty()) {
			assert(try_catch.HasCaught());
			k8_exception(&try_catch);
			return false;
		} else {
			assert(!try_catch.HasCaught());
			if (prt_rst && !result->IsUndefined()) {
				v8::String::AsciiValue str(result);
				puts(k8_cstr(str));
			}
			return true;
		}
	}
}

v8::Handle<v8::String> k8_readfile(const char *name) // Read the entire file. Copied from v8/shell.cc
{
	FILE* file = fopen(name, "rb");
	if (file == NULL) {
		fprintf(stderr, "ERROR: fail to open file '%s'.\n", name);
		return v8::Handle<v8::String>();
	}

	fseek(file, 0, SEEK_END);
	int size = ftell(file);
	rewind(file);

	char* chars = new char[size + 1];
	chars[size] = '\0';
	for (int i = 0; i < size;) {
		int read = static_cast<int>(fread(&chars[i], 1, size - i, file));
		i += read;
	}
	fclose(file);

	if (size > 2 && strncmp(chars, "#!", 2) == 0) { // then skip the "#!" line
		int i;
		for (i = 0; i < size; ++i)
			if (chars[i] == '\n') break;
		size -= i + 1;
		memmove(chars, &chars[i+1], size);
	}

	v8::Handle<v8::String> result = v8::String::New(chars, size);
	delete[] chars;
	return result;
}

/******************************
 *** New built-in functions ***
 ******************************/

JS_METHOD(k8_print, args) // print(): print to stdout; TAB demilited if multiple arguments are provided
{
	for (int i = 0; i < args.Length(); i++) {
		v8::HandleScope handle_scope;
		if (i) putchar('\t');
		v8::String::AsciiValue str(args[i]);
		fputs(k8_cstr(str), stdout);
	}
	putchar('\n');
	return v8::Undefined();
}

JS_METHOD(k8_warn, args) // print(): print to stdout; TAB demilited if multiple arguments are provided
{
	for (int i = 0; i < args.Length(); i++) {
		v8::HandleScope handle_scope;
		if (i) fputc('\t', stderr);
		v8::String::AsciiValue str(args[i]);
		fputs(k8_cstr(str), stderr);
	}
	fputc('\n', stderr);
	return v8::Undefined();
}

JS_METHOD(k8_exit, args) // exit()
{
	int exit_code = args[0]->Int32Value();
	fflush(stdout); fflush(stderr);
	exit(exit_code);
	return v8::Undefined();
}

JS_METHOD(k8_load, args) // load(): Load and execute a JS file. It also searches ONE path in $K8_LIBRARY_PATH
{
	char buf[1024], *path = getenv("K8_PATH");
	FILE *fp;
	for (int i = 0; i < args.Length(); i++) {
		v8::HandleScope handle_scope;
		v8::String::Utf8Value file(args[i]);
		buf[0] = 0;
		if ((fp = fopen(*file, "r")) != 0) {
			fclose(fp);
			strcpy(buf, *file);
		} else if (path) { // TODO: to allow multiple paths separated by ":"
			strcpy(buf, path); strcat(buf, "/"); strcat(buf, *file);
			if ((fp = fopen(buf, "r")) == 0) buf[0] = 0;
			else fclose(fp);
		}
		if (buf[0] == 0) return JS_THROW(Error, "[load] fail to locate the file");
		v8::Handle<v8::String> source = k8_readfile(buf);
		if (!k8_execute(source, v8::String::New(*file), false))
			return JS_THROW(Error, "[load] fail to execute the file");
	}
	return v8::Undefined();
}

/********************
 *** Bytes object ***
 ********************/

typedef struct {
	int32_t n, m, tshift;
	v8::ExternalArrayType eta;
	uint8_t *a;
} kvec8_t;

static inline void set_length(const v8::Handle<v8::Object> &obj, const kvec8_t *v)
{
	obj->SetIndexedPropertiesToExternalArrayData(v->a, v->eta, v->n >> v->tshift);
}

static inline void kv_set_type(kvec8_t *v, const char *type)
{
	if (type == 0) v->tshift = 0, v->eta = v8::kExternalUnsignedByteArray;
	else if (strcmp(type, "int8_t") == 0) v->tshift = 0, v->eta = v8::kExternalByteArray;
	else if (strcmp(type, "uint8_t") == 0) v->tshift = 0, v->eta = v8::kExternalUnsignedByteArray;
	else if (strcmp(type, "int16_t") == 0) v->tshift = 1, v->eta = v8::kExternalShortArray;
	else if (strcmp(type, "uint16_t") == 0) v->tshift = 1, v->eta = v8::kExternalUnsignedShortArray;
	else if (strcmp(type, "int32_t") == 0) v->tshift = 2, v->eta = v8::kExternalIntArray;
	else if (strcmp(type, "uint32_t") == 0) v->tshift = 2, v->eta = v8::kExternalUnsignedIntArray;
	else if (strcmp(type, "float") == 0) v->tshift = 2, v->eta = v8::kExternalFloatArray;
	else if (strcmp(type, "double") == 0) v->tshift = 3, v->eta = v8::kExternalDoubleArray;
	else v->tshift = 0, v->eta = v8::kExternalUnsignedByteArray;
}

JS_METHOD(k8_bytes, args)
{
	v8::HandleScope scope;
	ASSERT_CONSTRUCTOR(args);
	kvec8_t *a;
	a = (kvec8_t*)calloc(1, sizeof(kvec8_t));
	kv_set_type(a, 0);
	if (args.Length() > 1 && args[1]->IsString()) {
		v8::String::AsciiValue type(args[1]);
		kv_set_type(a, *type);
	}
	if (args.Length()) {
		a->m = a->n = args[0]->Int32Value() << a->tshift;
		a->a = (uint8_t*)calloc(a->n, 1); // NB: we are expecting malloc/calloc/realloc() only allocate aligned memory
	}
	set_length(args.This(), a);
	SAVE_PTR(args, 0, a);
	return args.This();
}

static inline void kv_recapacity(kvec8_t *a, int32_t m)
{
	kroundup32(m);
	if (a->m != m) {
		a->a = (uint8_t*)realloc(a->a, m);
		if (a->m < m) memset(&a->a[a->m], 0, m - a->m);
		v8::V8::AdjustAmountOfExternalAllocatedMemory(m - a->m);
		a->m = m;
	}
	if (a->n > a->m) a->n = a->m;
}

JS_METHOD(k8_bytes_cast, args)
{
	v8::HandleScope scope;
	kvec8_t *v = LOAD_PTR(args, 0, kvec8_t*);
	if (args.Length()) {
		v8::String::AsciiValue type(args[0]);
		kv_set_type(v, *type);
	} else kv_set_type(v, 0);
	set_length(args.This(), v);
	return v8::Undefined();
}

JS_METHOD(k8_bytes_destroy, args)
{
	v8::HandleScope scope;
	kvec8_t *a = LOAD_PTR(args, 0, kvec8_t*);
	free(a->a);
	v8::V8::AdjustAmountOfExternalAllocatedMemory(-a->m);
	a->a = 0; a->n = a->m = 0;
	set_length(args.This(), a);
	free(a);
	SAVE_PTR(args, 0, 0);
	return v8::Undefined();
}

JS_METHOD(k8_bytes_set, args)
{
#define _extend_vec_(_l_) do { \
		if (pos + (int32_t)(_l_) >= a->m) \
			kv_recapacity(a, pos + (_l_)); \
		if (pos + (int32_t)(_l_) >= a->n >> a->tshift << a->tshift) { \
			a->n = pos + (_l_); \
			set_length(args.This(), a); \
		} \
		cnt = (_l_); \
	} while (0)

	v8::HandleScope scope;
	kvec8_t *a = LOAD_PTR(args, 0, kvec8_t*);
	if (args.Length() == 0) return v8::Undefined();
	int cnt = 0;
	int32_t pos = args.Length() >= 2? args[1]->Int32Value() : a->n>>a->tshift;
	if (args[0]->IsNumber()) {
		_extend_vec_(1<<a->tshift);
		if (a->eta == v8::kExternalUnsignedByteArray) a->a[pos] = args[0]->Uint32Value();
		else if (a->eta == v8::kExternalDoubleArray) ((double*)a->a)[pos] = args[0]->NumberValue();
		else if (a->eta == v8::kExternalFloatArray) ((float*)a->a)[pos] = args[0]->NumberValue();
		else if (a->eta == v8::kExternalByteArray) ((int8_t*)a->a)[pos] = args[0]->Int32Value();
		else if (a->eta == v8::kExternalIntArray) ((int32_t*)a->a)[pos] = args[0]->Int32Value();
		else if (a->eta == v8::kExternalUnsignedIntArray) ((uint32_t*)a->a)[pos] = args[0]->Uint32Value();
		else if (a->eta == v8::kExternalShortArray) ((int16_t*)a->a)[pos] = args[0]->Int32Value();
		else if (a->eta == v8::kExternalUnsignedShortArray) ((uint16_t*)a->a)[pos] = args[0]->Uint32Value();
	} else if (args[0]->IsString()) {
		v8::String::AsciiValue str(args[0]);
		const char *cstr = *str;
		_extend_vec_(str.length());
		for (int i = 0; i < str.length(); ++i) a->a[i+pos] = uint8_t(cstr[i]);
	} else if (args[0]->IsArray()) {
		unsigned i;
		v8::Handle<v8::Array> array = v8::Handle<v8::Array>::Cast(args[0]);
		_extend_vec_(array->Length()<<a->tshift);
		if (a->eta == v8::kExternalUnsignedByteArray) for (i = 0; i < array->Length(); ++i) a->a[pos + i] = array->Get(v8::Integer::New(i))->Uint32Value();
		else if (a->eta == v8::kExternalDoubleArray) for (i = 0; i < array->Length(); ++i) ((double*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->NumberValue();
		else if (a->eta == v8::kExternalFloatArray) for (i = 0; i < array->Length(); ++i) ((float*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->NumberValue();
		else if (a->eta == v8::kExternalByteArray) for (i = 0; i < array->Length(); ++i) ((int8_t*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->Int32Value();
		else if (a->eta == v8::kExternalIntArray) for (i = 0; i < array->Length(); ++i) ((int32_t*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->Int32Value();
		else if (a->eta == v8::kExternalUnsignedIntArray) for (i = 0; i < array->Length(); ++i) ((uint32_t*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->Uint32Value();
		else if (a->eta == v8::kExternalShortArray) for (i = 0; i < array->Length(); ++i) ((int16_t*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->Int32Value();
		else if (a->eta == v8::kExternalUnsignedShortArray) for (i = 0; i < array->Length(); ++i) ((uint16_t*)a->a)[pos + i] = array->Get(v8::Integer::New(i))->Uint32Value();
	} else if (args[0]->IsObject()) {
		v8::Handle<v8::Object> b = v8::Handle<v8::Object>::Cast(args[0]); // TODO: check b is a 'Bytes' instance
		kvec8_t *a2 = reinterpret_cast<kvec8_t*>(b->GetAlignedPointerFromInternalField(0));
		_extend_vec_(a2->n);
		memcpy(a->a + (pos << a->tshift), a2->a, a2->n);
	}
	return scope.Close(v8::Integer::New(cnt));
}

JS_METHOD(k8_bytes_toString, args)
{
	v8::HandleScope scope;
	kvec8_t *a = LOAD_PTR(args, 0, kvec8_t*);
	return scope.Close(v8::String::New((char*)a->a, a->n));
}

v8::Handle<v8::Value> k8_bytes_length_getter(v8::Local<v8::String> property, const v8::AccessorInfo &info)
{
	kvec8_t *a = LOAD_PTR(info, 0, kvec8_t*);
	return v8::Integer::New(a->n >> a->tshift);
}

void k8_bytes_length_setter(v8::Local<v8::String> property, v8::Local<v8::Value> value, const v8::AccessorInfo &info)
{
	kvec8_t *a = LOAD_PTR(info, 0, kvec8_t*);
	int32_t n_old = a->n;
	a->n = value->Int32Value() << a->tshift;
	if (a->n > a->m) kv_recapacity(a, a->n);
	if (n_old != a->n) set_length(info.This(), a);
}

v8::Handle<v8::Value> k8_bytes_capacity_getter(v8::Local<v8::String> property, const v8::AccessorInfo &info)
{
	kvec8_t *a = LOAD_PTR(info, 0, kvec8_t*);
	return v8::Integer::New(a->m >> a->tshift);
}

void k8_bytes_capacity_setter(v8::Local<v8::String> property, v8::Local<v8::Value> value, const v8::AccessorInfo &info)
{
	kvec8_t *a = LOAD_PTR(info, 0, kvec8_t*);
	int32_t n_old = a->n;
	kv_recapacity(a, value->Int32Value() << a->tshift);
	if (n_old != a->n) set_length(info.This(), a);
}

/**********************************************
 *** Generic stream buffer from klib/kseq.h ***
 **********************************************/

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: " \n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

typedef struct __kstream_t {
	unsigned char *buf;
	int begin, end, is_eof, buf_size;
} kstream_t;

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

static inline kstream_t *ks_init(int __bufsize)
{
	kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));
	ks->buf_size = __bufsize;
	ks->buf = (unsigned char*)malloc(__bufsize);
	return ks;
}

static inline void ks_destroy(kstream_t *ks)
{ 
	if (ks) { free(ks->buf); free(ks); }
}

template<typename file_t, typename reader_t>
static int ks_getc(file_t &fp, kstream_t *ks, reader_t reader)
{
	if (ks->is_eof && ks->begin >= ks->end) return -1;
	if (ks->begin >= ks->end) {
		ks->begin = 0;
		ks->end = reader(fp, ks->buf, ks->buf_size);
		if (ks->end < ks->buf_size) ks->is_eof = 1;
		if (ks->end == 0) return -1;
	}
	return (int)ks->buf[ks->begin++];
}

template<typename file_t, typename reader_t>
static size_t ks_read(file_t &fp, kstream_t *ks, uint8_t *buf, long len, reader_t reader)
{
	long offset = 0;
	if (ks->is_eof && ks->begin >= ks->end) return -1;
	while (len > ks->end - ks->begin) {
		int l = ks->end - ks->begin;
		memcpy(buf + offset, ks->buf + ks->begin, l);
		len -= l; offset += l;
		ks->begin = 0;
		ks->end = reader(fp, ks->buf, ks->buf_size);
		if (ks->end < ks->buf_size) ks->is_eof = 1;
		if (ks->end == 0) return offset;
	}
	memcpy(buf + offset, ks->buf + ks->begin, len);
	ks->begin += len;
	return offset + len;
}

template<typename file_t, typename reader_t>
static int ks_getuntil(file_t &fp, kstream_t *ks, kvec8_t *kv, int delimiter, int *dret, int offset, reader_t reader)
{
	if (dret) *dret = 0;
	kv->n = offset >= 0? offset : 0;
	if (ks->begin >= ks->end && ks->is_eof) return -1;
	for (;;) {
		int i;
		if (ks->begin >= ks->end) {
			if (!ks->is_eof) {
				ks->begin = 0;
				ks->end = reader(fp, ks->buf, ks->buf_size);
				if (ks->end < ks->buf_size) ks->is_eof = 1;
				if (ks->end == 0) break;
			} else break;
		}
		if (delimiter == KS_SEP_LINE) {
			for (i = ks->begin; i < ks->end; ++i)
				if (ks->buf[i] == '\n') break;
		} else if (delimiter > KS_SEP_MAX) {
			for (i = ks->begin; i < ks->end; ++i)
				if (ks->buf[i] == delimiter) break;
		} else if (delimiter == KS_SEP_SPACE) {
			for (i = ks->begin; i < ks->end; ++i)
				if (isspace(ks->buf[i])) break;
		} else if (delimiter == KS_SEP_TAB) {
			for (i = ks->begin; i < ks->end; ++i)
				if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break;
		} else i = 0; /* never come to here! */
		if (kv->m - kv->n < i - ks->begin + 1) {
			kv->m = kv->n + (i - ks->begin) + 1;
			kroundup32(kv->m);
			kv->a = (uint8_t*)realloc(kv->a, kv->m);
		}
		memcpy(kv->a + kv->n, ks->buf + ks->begin, i - ks->begin);
		kv->n = kv->n + (i - ks->begin);
		ks->begin = i + 1;
		if (i < ks->end) {
			if (dret) *dret = ks->buf[i];
			break;
		}
	}
	if (kv->a == 0) {
		kv->m = 1;
		kv->a = (uint8_t*)calloc(1, 1);
	} else if (delimiter == KS_SEP_LINE && kv->n > 1 && kv->a[kv->n-1] == '\r') --kv->n;
	kv->a[kv->n] = '\0';
	return kv->n;
}

#define KS_BUF_SIZE 0x10000

/*******************
 *** File object ***
 *******************/

JS_METHOD(k8_file, args) // new File(fileName=stdin, mode="r").
{
	v8::HandleScope scope;
	ASSERT_CONSTRUCTOR(args);
	int fd = -1;
	FILE *fpw = 0;  // write ordinary files
	gzFile fpr = 0; // zlib transparently reads both ordinary and zlib/gzip files
	if (args.Length()) {
		SAVE_VALUE(args, 0, args[0]); // InternalField[0] keeps the file name
		v8::String::AsciiValue file(args[0]);
		if (args[0]->IsUint32()) fd = args[0]->Int32Value();
		if (args.Length() >= 2) {
			SAVE_VALUE(args, 1, args[1]); // InternalField[1] keeps the mode
			v8::String::AsciiValue mode(args[1]);
			const char *cstr = k8_cstr(mode);
			if (strchr(cstr, 'w')) fpw = fd >= 0? fdopen(fd, cstr) : fopen(k8_cstr(file), cstr);
			else fpr = fd >= 0? gzdopen(fd, cstr) : gzopen(k8_cstr(file), cstr);
		} else {
			SAVE_VALUE(args, 1, JS_STR("r"));
			fpr = fd >= 0? gzdopen(fd, "r") : gzopen(*file, "r");
		}
		if (fpr == 0 && fpw == 0)
			return JS_THROW(Error, "[File] Fail to open the file");
	} else {
		SAVE_VALUE(args, 0, JS_STR("-"));
		SAVE_VALUE(args, 1, JS_STR("r"));
		fpr = gzdopen(fileno(stdin), "r");
	}
	SAVE_PTR(args, 2, fpr); // InternalField[2] keeps the reading file pointer
	SAVE_PTR(args, 3, fpw); // InternalField[3] keeps the writing file pointer
	if (fpr) {
		kstream_t *ks = ks_init(KS_BUF_SIZE);
		v8::V8::AdjustAmountOfExternalAllocatedMemory(KS_BUF_SIZE);
		SAVE_PTR(args, 4, ks);
	} else SAVE_PTR(args, 4, 0);
	return args.This();
}

JS_METHOD(k8_file_close, args) // File::close(). Close the file.
{
	gzFile fpr = LOAD_PTR(args, 2, gzFile);
	FILE  *fpw = LOAD_PTR(args, 3, FILE*);
	if (fpr) {
		gzclose(fpr);
		kstream_t *ks = LOAD_PTR(args, 4, kstream_t*);
		ks_destroy(ks);
		v8::V8::AdjustAmountOfExternalAllocatedMemory(-KS_BUF_SIZE);
	}
	if (fpw) fclose(fpw);
	SAVE_PTR(args, 2, 0); SAVE_PTR(args, 3, 0); SAVE_PTR(args, 4, 0);
	return v8::Undefined();
}

JS_METHOD(k8_file_read, args) // File::read(), read(buf, offset, length)
{
	v8::HandleScope scope;
	gzFile fp = LOAD_PTR(args, 2, gzFile);
	kstream_t *ks = LOAD_PTR(args, 4, kstream_t*);
	if (fp == 0) return JS_ERROR("file is not open for read");
	if (args.Length() == 0) { // read()
		int c = ks_getc(fp, ks, gzread);
		return scope.Close(v8::Integer::New(c));
	} else if (args.Length() == 3 && args[0]->IsObject() && args[1]->IsUint32() && args[2]->IsUint32()) { // read(buf, offset, length)
		long off = args[1]->Int32Value(), len = args[2]->Int32Value();
		v8::Handle<v8::Object> b = v8::Handle<v8::Object>::Cast(args[0]); // TODO: check b is a 'Bytes' instance
		kvec8_t *kv = reinterpret_cast<kvec8_t*>(b->GetAlignedPointerFromInternalField(0));
		if (len + off > kv->m) kv_recapacity(kv, len + off);
		len = ks_read(fp, ks, kv->a + off, len, gzread);
		if (len + off > kv->n) {
			kv->n = len + off;
			set_length(b, kv);
		}
		return scope.Close(v8::Integer::New(len));
	}
	return JS_ERROR("misused File.prototype.read()");
}

JS_METHOD(k8_file_write, args) // File::write(str). Write $str and return the number of written characters
{
	v8::HandleScope scope;
	FILE *fp = LOAD_PTR(args, 3, FILE*);
	if (fp == 0) return JS_ERROR("file is not open for write");
	if (args.Length() == 0) return scope.Close(v8::Integer::New(0));
	long len = 0;
	if (args[0]->IsString()) {
		v8::String::AsciiValue vbuf(args[0]);
		len = fwrite(*vbuf, 1, vbuf.length(), fp);
	} else if (args[0]->IsObject()) {
		v8::Handle<v8::Object> b = v8::Handle<v8::Object>::Cast(args[0]); // TODO: check b is a 'Bytes' instance
		kvec8_t *kv = reinterpret_cast<kvec8_t*>(b->GetAlignedPointerFromInternalField(0));
		len = fwrite(kv->a, 1, kv->n, fp);
	}
	return scope.Close(v8::Integer::New(len));
}

JS_METHOD(k8_file_readline, args) // see iStream::readline(sep=line) for details
{
	v8::HandleScope scope;
	gzFile fpr = LOAD_PTR(args, 2, gzFile);
	kstream_t *ks = LOAD_PTR(args, 4, kstream_t*);
	if (fpr == 0) return JS_ERROR("file is not open for read");
	if (!args.Length() || !args[0]->IsObject()) return v8::Null(); // TODO: when there are no parameters, skip a line
	v8::Handle<v8::Object> b = v8::Handle<v8::Object>::Cast(args[0]); // TODO: check b is a 'Bytes' instance
	kvec8_t *kv = reinterpret_cast<kvec8_t*>(b->GetAlignedPointerFromInternalField(0));
	int dret, ret, sep = KS_SEP_LINE;
	if (args.Length() > 1) { // by default, the delimitor is new line
		if (args[1]->IsString()) { // if 1st argument is a string, set the delimitor to the 1st charactor of the string
			v8::String::AsciiValue str(args[1]);
			sep = int(k8_cstr(str)[0]);
		} else if (args[1]->IsInt32()) // if 1st argument is an integer, set the delimitor to the integer: 0=>isspace(); 1=>isspace()&&!' '; 2=>newline
			sep = args[1]->Int32Value();
	}
	int offset = 0;
	if (args.Length() > 2) {
		if (args[2]->IsUint32()) offset = args[2]->Uint32Value();
		else if (args[2]->IsBoolean()) offset = args[2]->BooleanValue()? kv->n : 0;
	}
	ret = ks_getuntil(fpr, ks, kv, sep, &dret, offset, gzread);
	set_length(b, kv);
	return ret >= 0? scope.Close(v8::Integer::New(dret)) : scope.Close(v8::Integer::New(ret));
}

/******************
 *** Set object ***
 ******************/

#include "khash.h"
KHASH_MAP_INIT_STR(str, kh_cstr_t)
typedef khash_t(str) *strset_t;

static const char *k8_empty_str = "";

JS_METHOD(k8_map, args)
{
	v8::HandleScope scope;
	ASSERT_CONSTRUCTOR(args);
	strset_t h = kh_init(str);
	SAVE_PTR(args, 0, h);
	return args.This();
}

JS_METHOD(k8_map_put, args)
{
	v8::HandleScope scope;
	strset_t h = LOAD_PTR(args, 0, strset_t);
	if (args.Length()) {
		v8::String::AsciiValue s(args[0]);
		const char *cstr = k8_cstr(s);
		int absent;
		khint_t k = kh_put(str, h, cstr, &absent);
		if (absent) {
			kh_key(h, k) = strdup(cstr);
			kh_val(h, k) = k8_empty_str;
		}
		if (args.Length() > 1) {
			v8::String::AsciiValue v(args[1]);
			if (kh_val(h, k) != k8_empty_str)
				free((char*)kh_val(h, k));
			kh_val(h, k) = strdup(k8_cstr(v));
		}
	} else return JS_ERROR("misused Map.prototype.put()");
	return v8::Undefined();
}

JS_METHOD(k8_map_get, args)
{
	v8::HandleScope scope;
	strset_t h = LOAD_PTR(args, 0, strset_t);
	if (args.Length()) {
		v8::String::AsciiValue s(args[0]);
		const char *cstr = k8_cstr(s);
		khint_t k = kh_get(str, h, cstr);
		return k == kh_end(h)? v8::Null() : scope.Close(v8::String::New(kh_val(h, k), strlen(kh_val(h, k))));
	} else return JS_ERROR("misused Map.prototype.get()");
}

JS_METHOD(k8_map_del, args)
{
	v8::HandleScope scope;
	strset_t h = LOAD_PTR(args, 0, strset_t);
	if (args.Length()) {
		v8::String::AsciiValue s(args[0]);
		const char *cstr = k8_cstr(s);
		khint_t k = kh_get(str, h, cstr);
		if (k < kh_end(h)) {
			free((char*)kh_key(h, k));
			if (kh_val(h, k) != k8_empty_str)
				free((char*)kh_val(h, k));
			kh_del(str, h, k);
		}
	} else return JS_ERROR("misused Map.prototype.del()");
	return v8::Undefined();
}

JS_METHOD(k8_map_destroy, args)
{
	v8::HandleScope scope;
	strset_t h = LOAD_PTR(args, 0, strset_t);
	khint_t k;
	for (k = 0; k != kh_end(h); ++k)
		if (kh_exist(h, k)) {
			free((char*)kh_key(h, k));
			if (kh_val(h, k) != k8_empty_str)
				free((char*)kh_val(h, k));
		}
	kh_destroy(str, h);
	SAVE_PTR(args, 0, 0);
	return v8::Undefined();
}

/***********************
 *** Getopt from BSD ***
 ***********************/

// Modified from getopt.c from BSD, 3-clause BSD license. Copyright (c) 1987-2002 The Regents of the University of California.
// We do not use system getopt() because it may parse "-v" in "k8 prog.js -v".

int opterr = 1, optind = 1, optopt, optreset;
char *optarg;

int getopt(int nargc, char * const *nargv, const char *ostr)
{
	static char *place = 0;
	const char *oli;
	if (optreset || !place || !*place) {
		optreset = 0;
		if (optind >= nargc || *(place = nargv[optind]) != '-') {
			place = 0;
			return -1;
		}
		if (place[1] && *++place == '-') {
			++optind, place = 0;
			return -1;
		}
	}
	if ((optopt = *place++) == ':' || !(oli = strchr(ostr, optopt))) {
		if (optopt == '-') return -1;
		if (!*place) ++optind;
		if (opterr && *ostr != ':') fprintf(stderr, "%s: illegal option -- %c\n", __FILE__, optopt);
		return '?';
	}
	if (*++oli != ':') {
		optarg = 0;
		if (!*place) ++optind;
	} else {
		if (*place) optarg = place;
		else if (nargc <= ++optind) {
			place = 0;
			if (*ostr == ':') return ':';
			if (opterr) fprintf(stderr, "%s: option requires an argument -- %c\n", __FILE__, optopt);
			return '?';
		} else optarg = nargv[optind];
		place = 0;
		++optind;
	}
	return optopt;
}

/*********************
 *** Main function ***
 *********************/

static v8::Persistent<v8::Context> CreateShellContext() // adapted from shell.cc
{
	v8::Handle<v8::ObjectTemplate> global = v8::ObjectTemplate::New();
	global->Set(JS_STR("print"), v8::FunctionTemplate::New(k8_print));
	global->Set(JS_STR("warn"), v8::FunctionTemplate::New(k8_warn));
	global->Set(JS_STR("exit"), v8::FunctionTemplate::New(k8_exit));
	global->Set(JS_STR("load"), v8::FunctionTemplate::New(k8_load));
	{ // add the 'Bytes' object
		v8::HandleScope scope;
		v8::Handle<v8::FunctionTemplate> ft = v8::FunctionTemplate::New(k8_bytes);
		ft->SetClassName(JS_STR("Bytes"));

		v8::Handle<v8::ObjectTemplate> ot = ft->InstanceTemplate();
		ot->SetInternalFieldCount(1);
		ot->SetAccessor(JS_STR("length"), k8_bytes_length_getter, k8_bytes_length_setter, v8::Handle<v8::Value>(), v8::DEFAULT, static_cast<v8::PropertyAttribute>(v8::DontDelete));
		ot->SetAccessor(JS_STR("capacity"), k8_bytes_capacity_getter, k8_bytes_capacity_setter, v8::Handle<v8::Value>(), v8::DEFAULT, static_cast<v8::PropertyAttribute>(v8::DontDelete));

		v8::Handle<v8::ObjectTemplate> pt = ft->PrototypeTemplate();
		pt->Set("cast", v8::FunctionTemplate::New(k8_bytes_cast));
		pt->Set("set", v8::FunctionTemplate::New(k8_bytes_set));
		pt->Set("toString", v8::FunctionTemplate::New(k8_bytes_toString));
		pt->Set("destroy", v8::FunctionTemplate::New(k8_bytes_destroy));

		global->Set("Bytes", ft);	
	}
	{ // add the 'File' object
		v8::HandleScope scope;
		v8::Handle<v8::FunctionTemplate> ft = v8::FunctionTemplate::New(k8_file);
		ft->SetClassName(JS_STR("File"));
		ft->InstanceTemplate()->SetInternalFieldCount(5); // (fn, mode, fpr, fpw)
		v8::Handle<v8::ObjectTemplate> pt = ft->PrototypeTemplate();
		pt->Set("read", v8::FunctionTemplate::New(k8_file_read));
		pt->Set("readline", v8::FunctionTemplate::New(k8_file_readline));
		pt->Set("write", v8::FunctionTemplate::New(k8_file_write));
		pt->Set("close", v8::FunctionTemplate::New(k8_file_close));
		pt->Set("destroy", v8::FunctionTemplate::New(k8_file_close));
		global->Set("File", ft);	
	}
	{ // add the 'Set' object
		v8::HandleScope scope;
		v8::Handle<v8::FunctionTemplate> ft = v8::FunctionTemplate::New(k8_map);
		ft->SetClassName(JS_STR("Map"));
		ft->InstanceTemplate()->SetInternalFieldCount(1);
		v8::Handle<v8::ObjectTemplate> pt = ft->PrototypeTemplate();
		pt->Set("put", v8::FunctionTemplate::New(k8_map_put));
		pt->Set("get", v8::FunctionTemplate::New(k8_map_get));
		pt->Set("del", v8::FunctionTemplate::New(k8_map_del));
		pt->Set("destroy", v8::FunctionTemplate::New(k8_map_destroy));
		global->Set("Map", ft);	
	}
	return v8::Context::New(NULL, global);
}

int main(int argc, char* argv[])
{
	v8::V8::SetFlagsFromCommandLine(&argc, argv, true);

	// set --max_old_space_size. We have to do it before CreateShellContext(), I guess.
	int c;
	char flag_buf[256], *size = 0;
	while ((c = getopt(argc, argv, "ve:E:M:")) >= 0)
		if (c == 'M') size = optarg;
	strcat(strcpy(flag_buf, "--max_old_space_size="), size? size : "16384");
	v8::V8::SetFlagsFromString(flag_buf, strlen(flag_buf));
	opterr = optind = 1;

	int ret = 0;
	if (argc == 1) {
		fprintf(stderr, "Usage: k8 [-v] [-e jsSrc] [-E jsSrc] [-M maxRSS] <src.js> [arguments]\n");
		return 1;
	}
	{
		v8::HandleScope scope;
		v8::Persistent<v8::Context> context = CreateShellContext();
		if (context.IsEmpty()) {
			fprintf(stderr, "Error creating context\n");
			return 1;
		}
		context->Enter();
		while ((c = getopt(argc, argv, "ve:E:M:")) >= 0) // parse k8 related command line options
			if (c == 'e' || c == 'E') {
				if (!k8_execute(JS_STR(optarg), JS_STR("CLI"), (c == 'E'))) { // note the difference between 'e' and 'E'
					ret = 1;
					break;
				}
			} else if (c == 'v') printf("V8: %s\nK8: %s\n", v8::V8::GetVersion(), K8_VERSION);
		if (!ret && optind != argc) {
			v8::HandleScope scope2;
			v8::Local<v8::Array> args = v8::Array::New(argc - optind - 1);
			for (int i = optind + 1; i < argc; ++i)
				args->Set(v8::Integer::New(i - optind - 1), JS_STR(argv[i]));
			context->Global()->Set(JS_STR("arguments"), args);
			if (!k8_execute(k8_readfile(argv[optind]), JS_STR(argv[optind]), false)) ret = 1;
		}
		context->Exit();
		context.Dispose();
	}
	v8::V8::Dispose();
	return ret;
}

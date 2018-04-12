FAQ
---

#### 1. What is K8?

K8 is a Javascript shell based on Google's [V8 Javascript engine][1]. It adds
the support of flexible byte arrays and file I/O. K8 is implemented in one C++
source file. The only dependency is zlib in addition to V8.

#### 2. There are many Javascript shells with much richer features. What makes K8 special?

To some extent, [Node.js][2], [Narwhal][3], [SilkJS][4], [TeaJS][5] and
[Sorrow.js][6] are all Javascript shells. They not only provide binary storage
and file I/O, the features available in K8, but also implement much richer
functionality such as network I/O and database binding. However, most of the
existing Javascript shells are designed for server-side applications, but not
for general use cases as we do with Perl/Ruby/Python.  Take the popular Node.js
as an example. Node.js mixes file I/O and file system operations, two distinct
concepts, in one [File System module][7].  In the module, we can either read an
entire file or a fixed-length data blob, but are unable to read a line as is
provided by most other programming languages. Many other JS shell
implementations follow the [CommonJS APIs][9], which have a similar problem: no
usable APIs for general-purpose file I/O. After all these efforts, on file I/O,
we even do not have a JS shell matching the usability of C, let alone
high-level programming languages such as Perl and Python.

K8 aims to provide C-like file I/O APIs. It adds a `File` object for buffered
file reading, a `Bytes` object for flexible binary storage and a `Map` object
for a hash map without hitting the memory limit of V8.

#### 3. How to compile K8? Are there compiled binaries?

You need to first compile V8 and then compile and link K8. Here is the full procedure:

```sh
# download compilable V8 source code; K8 only works with v8-3.16
wget -O- https://github.com/attractivechaos/k8/releases/download/v0.2.1/v8-3.16.4.tar.bz2 | tar jxf -
# compile V8
cd v8-3.16.4 && make -j4 x64.release
# compile K8
g++ -O2 -Wall -o k8 -Iinclude ../k8.cc -lpthread -lz `find out -name "libv8_base.a"` `find out -name "libv8_snapshot.a"`
```

Alternatively, you may download the precompiled binaries for Mac and Linux from
the [release page][release].

#### 4. An earlier version of K8 implemented a generic buffered stream. Why has it been removed?

To implement a generic buffered stream, we need to call a Javascript `read`
function in C++ and transform between Javascript and C++ data representation.
This procedure adds significant overhead. For the best performance on file
I/O, all the `iStream` functionality has been moved to `File`. Anyway, it
is not hard to implement buffered stream purely in Javascript.


API Documentations
------------------

All the following objects manage some memory outside the V8 garbage collector.
It is important to call the `close()` or the `destroy()` methods to deallocate
the memory to avoid memory leaks.

### Example

```javascript
var x = new Bytes(), y = new Bytes();
x.set('foo'); x.set([0x20,0x20]); x.set('bar'); x.set('F', 0); x[3]=0x2c;
print(x.toString())   // output: 'Foo, bar'
y.set('BAR'); x.set(y, 5)
print(x)              // output: 'Foo, BAR'
x.destroy(); y.destroy()
if (arguments.length) { // read and print file
  var x = new Bytes(), s = new File(arguments[0]);
  while (s.readline(x) >= 0) print(x)
  s.close(); x.destroy();
}
```

### The Bytes Object

`Bytes` provides a byte array. It has the following methods:

```javascript
// Create an array of type $type in length $len. $type can be: int8_t, uint8_t, int16_t,
// uint16_t, int32_t, uint32_t, float or double.
new Bytes(len, type)

// Equivalent to 'new Bytes(len, "uint8_t")'
new Bytes(len)

// Equivalent to 'new Bytes(0, "uint8_t")'
new Bytes()

// Property: get/set length of the array
.length

// Property: get/set the max capacity of the array
.capacity

// The index operator. If $pos goes beyond .length, undefined will be returned.
int obj[pos]

// Change the array type to $type, equivalent to changing the pointer type. .length and
// .capacity may be changed if the size of element is changed.
Bytes.prototype.cast(type)

// Equivalent to 'Bytes.prototype.cast("uint8_t")'
Bytes.prototype.cast()

// Deallocate the array. This is necessary as the memory is not managed by the V8 GC.
Bytes.prototype.destroy()

// Replace the byte array starting from $offset to $data, where $data can be a number,
// a string, an array or Bytes. The size of the array is modified if the new array
// is larger. Return the number of modified bytes. If only one byte needs to be
// changed, using the [] operator gives better performance.
int Bytes.prototype.set(data, offset)

// Append $data to the byte array
int Bytes.prototype.set(data)

// Convert the byte array to string
Bytes.prototype.toString()
```

### The File Object

`File` provides buffered file I/O. It has the following methods:

```javascript
// Open $fileName under $mode. $mode is in the same syntax as fopen(). Integer $fileName for
// a file descriptor. In particular, 0 for STDIN, 1 for STDOUT and 2 for STDERR.
new File(fileName, mode)

// Equivalent to 'new File(fileName, "r")'
new File(fileName)

// Equivalent to 'new File(0)'
new File()

// Read a byte. Return -1 if reaching end-of-file
int File.prototype.read()

// Read maximum $len bytes of data to $buf, starting from $offset. Return the number of
// bytes read to $buf. The size of $buf is unchanged unless it is smaller than $offset+$len.
int File.prototype.read(buf, offset, len)

// Write $data, which can be a string or Bytes(). Return the number of written bytes.
// This method replies on C's fwrite() for buffering.
int File.prototype.write(data)

// Read a line to $bytes starting from $offset, using $sep as the separator. $sep==0 sets
// the separator to isspace(), $sep==1 to (isspace() && !' ') and $sep==2 to newline. If
// $sep is a string, the first character in the string is the separator. Return the line
// length or -1 if reaching end-of-file.
int File.prototype.readline(bytes, sep, offset)

// Equivalent to 'File.prototype.readline(bytes, sep, 0)'
int File.prototype.readline(bytes, sep)

// Equivalent to 'File.prototype.readline(bytes, 2, 0)'
int File.prototype.readline(bytes)

// Close the file
File.prototype.close()
```

### The Map Object

`Map` provides a hash map implementation without using memory managed by V8. This can be helpful
when we want to stage a huge hash table in memory. `Map` has the following methods:

```javascript
// Initialize a hash map
new Map()

// Put a key-value string pair to a hash map
Map.prototype.put(key, value)

// Equivalent to 'Map.prototype.put(key, "")'
Map.prototype.put(key)

// Get a key. Return 'null' if 'key' is non-existing
string Map.prototype.get(key)

// Delete a key.
Map.prototype.del(key)

// Deallocate memory
Map.prototype.destroy()
```

[1]: http://code.google.com/p/v8/
[2]: http://nodejs.org/
[3]: https://github.com/tlrobinson/narwhal
[4]: http://silkjs.net/
[5]: http://code.google.com/p/teajs/
[6]: https://github.com/samlecuyer/sorrow.js
[7]: http://nodejs.org/api/fs.html
[8]: http://nodejs.org/api/stream.html
[9]: http://www.commonjs.org/specs/
[11]: https://sourceforge.net/projects/lh3/files/
[gyp]: https://gyp.gsrc.io/
[release]: https://github.com/attractivechaos/k8/releases

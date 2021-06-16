ref_data=$1
preset=$2

make clean && make no_opt=1
touch temp_read.fastq
./minimap2  -ax $2 $1 temp_read.fastq -Z 1 >/dev/null      


kv_file=$1"_"$2"_minimizers_key_value_sorted"  

full_path=`readlink -f $kv_file`

cd ./ext/TAL
make lisa_hash
./build-lisa-hash-index $full_path

rm ../../temp_read.fastq

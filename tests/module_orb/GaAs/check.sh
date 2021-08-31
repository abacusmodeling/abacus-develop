cd Table_SR0
for i in *dat
do
	diff $i ../ref/Table_SR0/$i
done
cd ../Table_TR0
for i in *dat
do
	diff $i ../ref/Table_TR0/$i
done
cd ..


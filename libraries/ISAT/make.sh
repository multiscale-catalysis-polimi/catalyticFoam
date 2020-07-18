mkdir lib

echo "Compiling ISAT..."

g++ -Iinclude -I$EIGEN_LIBRARY_PATH  -g -fpic -c include/ISAT.C include/binaryTree.C include/binaryNode.C include/chemComp.C

g++ -shared ISAT.o binaryTree.o binaryNode.o chemComp.o -o lib/libISAT.so 

ar rcs lib/libISAT.a ISAT.o binaryTree.o binaryNode.o chemComp.o

rm -f ISAT.o binaryTree.o binaryNode.o chemComp.o

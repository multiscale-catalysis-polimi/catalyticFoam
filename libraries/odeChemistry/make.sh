mkdir lib

export BOOST_LIBS='-Wl,--start-group -Wl,-Bstatic -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -Wl,--end-group -Wl,-Bdynamic'

echo "Compiling OdeChemistryNoMKL..."
g++ -w -g -Iinclude -I../OpenSMOKEpp4catalyticFoam/include -I$EIGEN_LIBRARY_PATH -I$BOOST_LIBRARY_PATH/include/ -I$RAPIDXML_LIBRARY_PATH -L$BOOST_LIBRARY_PATH/lib/ -lgfortran -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS $MKL_SUPPORT -fpic -c include/BatchReactorHeterogeneousConstantPressure.C include/BatchReactorHeterogeneousConstantVolume.C include/BatchReactorHomogeneousConstantPressure.C include/BatchReactorHomogeneousConstantVolume.C


g++ -shared BatchReactorHeterogeneousConstantPressure.o BatchReactorHeterogeneousConstantVolume.o BatchReactorHomogeneousConstantPressure.o BatchReactorHomogeneousConstantVolume.o -o lib/libOdeChemistryNoMKL.so 

ar rcs lib/libOdeChemistryNoMKL.a BatchReactorHeterogeneousConstantPressure.o BatchReactorHeterogeneousConstantVolume.o BatchReactorHomogeneousConstantPressure.o BatchReactorHomogeneousConstantVolume.o

rm -f BatchReactorHeterogeneousConstantPressure.o BatchReactorHeterogeneousConstantVolume.o BatchReactorHomogeneousConstantPressure.o BatchReactorHomogeneousConstantVolume.o

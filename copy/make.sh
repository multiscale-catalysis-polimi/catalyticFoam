g++ -w -static-libgcc -static-libstdc++ -std=c++11  -I/home/mauro/Simulations/Solver/catalyticPimpleFoam09_newOS/of5x/catalyticPimpleFoam/libraries/OpenSMOKEpp/source -I$EIGEN_LIBRARY_PATH -I$BOOST_LIBRARY_PATH/include/ -I$RAPIDXML_LIBRARY_PATH -o OpenSMOKE_CHEMKIN_preprocessor CHEMKIN_PreProcessor.C -L$BOOST_LIBRARY_PATH/lib/  -Wl,-Bstatic -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -Wl,-Bdynamic 



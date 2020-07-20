mkdir lib

export BOOST_LIBS='-Wl,--start-group -Wl,-Bstatic -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -Wl,--end-group -Wl,-Bdynamic'

echo "Compiling OpenSMOKEppDictionaries..."
g++ -w -Iinclude -I$EIGEN_LIBRARY_PATH -I$BOOST_LIBRARY_PATH/include/ -I$RAPIDXML_LIBRARY_PATH -L$BOOST_LIBRARY_PATH/lib/ -lgfortran -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS $MKL_SUPPORT -fpic -c include/dictionary/InputFileDictionary.C include/dictionary/OpenSMOKE_Dictionary.C include/dictionary/OpenSMOKE_DictionaryFile.C include/dictionary/OpenSMOKE_DictionaryGrammar.C include/dictionary/OpenSMOKE_DictionaryKeyWord.C include/dictionary/OpenSMOKE_DictionaryManager.C

g++ -shared InputFileDictionary.o OpenSMOKE_Dictionary.o OpenSMOKE_DictionaryFile.o OpenSMOKE_DictionaryGrammar.o OpenSMOKE_DictionaryKeyWord.o OpenSMOKE_DictionaryManager.o -o lib/libOpenSMOKEppDictionaries.so 

ar rcs lib/libOpenSMOKEppDictionaries.a InputFileDictionary.o OpenSMOKE_Dictionary.o OpenSMOKE_DictionaryFile.o OpenSMOKE_DictionaryGrammar.o OpenSMOKE_DictionaryKeyWord.o OpenSMOKE_DictionaryManager.o

rm -f InputFileDictionary.o OpenSMOKE_Dictionary.o OpenSMOKE_DictionaryFile.o OpenSMOKE_DictionaryGrammar.o OpenSMOKE_DictionaryKeyWord.o OpenSMOKE_DictionaryManager.o

echo "Compiling OpenSMOKEppKernel..."
g++ -w -Iinclude -I$EIGEN_LIBRARY_PATH -I$BOOST_LIBRARY_PATH/include/ -I$RAPIDXML_LIBRARY_PATH -L$BOOST_LIBRARY_PATH/lib/ -lgfortran -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS $MKL_SUPPORT  -fpic -c include/kernel/kinetics/ChebyshevPolynomialRateExpression.C include/kernel/kinetics/ExtendedFallOff.C include/kernel/kinetics/ExtendedPressureLogarithmicRateExpression.C include/kernel/kinetics/KineticsUtilityFunctions.C include/kernel/kinetics/PressureLogarithmicRateExpression.C include/kernel/kinetics/ReactionPolicy_CHEMKIN.C include/kernel/kinetics/ReactionPolicy_Solid_CHEMKIN.C include/kernel/kinetics/ReactionPolicy_Surface_CHEMKIN.C include/kernel/thermo/AtomicComposition.C include/kernel/thermo/AtomicElement.C include/kernel/thermo/AtomicElementMap.C include/kernel/thermo/InputFileCHEMKIN.C include/kernel/thermo/Species.C include/kernel/thermo/ThermoPolicy_CHEMKIN.C include/kernel/thermo/ThermoReader.C include/kernel/thermo/ThermoReaderPolicy_CHEMKIN.C include/kernel/transport/TransportPolicy_CHEMKIN.C include/kernel/transport/TransportReader.C include/kernel/transport/TransportReaderPolicy_CHEMKIN.C

g++ -shared ChebyshevPolynomialRateExpression.o ExtendedFallOff.o ExtendedPressureLogarithmicRateExpression.o KineticsUtilityFunctions.o PressureLogarithmicRateExpression.o ReactionPolicy_CHEMKIN.o ReactionPolicy_Solid_CHEMKIN.o ReactionPolicy_Surface_CHEMKIN.o AtomicComposition.o AtomicElement.o AtomicElementMap.o InputFileCHEMKIN.o Species.o ThermoPolicy_CHEMKIN.o ThermoReader.o ThermoReaderPolicy_CHEMKIN.o TransportPolicy_CHEMKIN.o TransportReader.o TransportReaderPolicy_CHEMKIN.o -o lib/libOpenSMOKEppKernel.so 

ar rcs lib/libOpenSMOKEppKernel.a ChebyshevPolynomialRateExpression.o ExtendedFallOff.o ExtendedPressureLogarithmicRateExpression.o KineticsUtilityFunctions.o PressureLogarithmicRateExpression.o ReactionPolicy_CHEMKIN.o ReactionPolicy_Solid_CHEMKIN.o ReactionPolicy_Surface_CHEMKIN.o AtomicComposition.o AtomicElement.o AtomicElementMap.o InputFileCHEMKIN.o Species.o ThermoPolicy_CHEMKIN.o ThermoReader.o ThermoReaderPolicy_CHEMKIN.o TransportPolicy_CHEMKIN.o TransportReader.o TransportReaderPolicy_CHEMKIN.o

rm -f ChebyshevPolynomialRateExpression.o ExtendedFallOff.o ExtendedPressureLogarithmicRateExpression.o KineticsUtilityFunctions.o PressureLogarithmicRateExpression.o ReactionPolicy_CHEMKIN.o ReactionPolicy_Solid_CHEMKIN.o ReactionPolicy_Surface_CHEMKIN.o AtomicComposition.o AtomicElement.o AtomicElementMap.o InputFileCHEMKIN.o Species.o ThermoPolicy_CHEMKIN.o ThermoReader.o ThermoReaderPolicy_CHEMKIN.o TransportPolicy_CHEMKIN.o TransportReader.o TransportReaderPolicy_CHEMKIN.o

echo "Compiling OpenSMOKEppMaps..."
g++ -w -Iinclude -I$EIGEN_LIBRARY_PATH -I$BOOST_LIBRARY_PATH/include/ -I$RAPIDXML_LIBRARY_PATH -L$BOOST_LIBRARY_PATH/lib/ -lgfortran -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS $MKL_SUPPORT  -fpic -c include/maps/FluxAnalysisMap.C include/maps/JacobianSparsityPatternMap.C include/maps/KineticsMap.C include/maps/KineticsMap_CHEMKIN.C include/maps/KineticsMap_Solid_CHEMKIN.C include/maps/KineticsMap_Surface_CHEMKIN.C include/maps/SensitivityMap.C include/maps/StoichiometricMap.C include/maps/ThermodynamicsMap.C include/maps/ThermodynamicsMap_CHEMKIN.C include/maps/ThermodynamicsMap_Solid_CHEMKIN.C include/maps/ThermodynamicsMap_Surface_CHEMKIN.C include/maps/TransportPropertiesMap.C include/maps/TransportPropertiesMap_CHEMKIN.C

g++ -shared FluxAnalysisMap.o JacobianSparsityPatternMap.o KineticsMap.o KineticsMap_CHEMKIN.o KineticsMap_Solid_CHEMKIN.o KineticsMap_Surface_CHEMKIN.o SensitivityMap.o StoichiometricMap.o ThermodynamicsMap.o ThermodynamicsMap_CHEMKIN.o ThermodynamicsMap_Solid_CHEMKIN.o ThermodynamicsMap_Surface_CHEMKIN.o TransportPropertiesMap.o TransportPropertiesMap_CHEMKIN.o -o lib/libOpenSMOKEppMaps.so 

ar rcs lib/libOpenSMOKEppMaps.a FluxAnalysisMap.o JacobianSparsityPatternMap.o KineticsMap.o KineticsMap_CHEMKIN.o KineticsMap_Solid_CHEMKIN.o KineticsMap_Surface_CHEMKIN.o SensitivityMap.o StoichiometricMap.o ThermodynamicsMap.o ThermodynamicsMap_CHEMKIN.o ThermodynamicsMap_Solid_CHEMKIN.o ThermodynamicsMap_Surface_CHEMKIN.o TransportPropertiesMap.o TransportPropertiesMap_CHEMKIN.o

rm -f FluxAnalysisMap.o JacobianSparsityPatternMap.o KineticsMap.o KineticsMap_CHEMKIN.o KineticsMap_Solid_CHEMKIN.o KineticsMap_Surface_CHEMKIN.o SensitivityMap.o StoichiometricMap.o ThermodynamicsMap.o ThermodynamicsMap_CHEMKIN.o ThermodynamicsMap_Solid_CHEMKIN.o ThermodynamicsMap_Surface_CHEMKIN.o TransportPropertiesMap.o TransportPropertiesMap_CHEMKIN.o

echo "Compiling OpenSMOKEppMath..."
g++ -w -Iinclude -fpermissive -I$EIGEN_LIBRARY_PATH -I$BOOST_LIBRARY_PATH/include/ -I$RAPIDXML_LIBRARY_PATH -L$BOOST_LIBRARY_PATH/lib/ -lgfortran -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS $MKL_SUPPORT  -fpic -c include/math/OpenSMOKEBaseClass.C include/math/OpenSMOKEClass.C include/math/OpenSMOKEStdInclude.C include/math/OpenSMOKEFunctions.C include/math/Conversions.C include/math/OpenSMOKEUtilities.C include/math/OpenSMOKEVector.C include/math/OpenSMOKEMatrix.C include/math/OpenSMOKEBandMatrix.C include/math/OpenSMOKE_MatrixSparsityPattern.C include/math/OpenSMOKELoad.C

g++ -shared OpenSMOKEBaseClass.o OpenSMOKEClass.o OpenSMOKEStdInclude.o OpenSMOKEFunctions.o Conversions.o OpenSMOKEUtilities.o OpenSMOKEVector.o OpenSMOKEMatrix.o OpenSMOKEBandMatrix.o OpenSMOKE_MatrixSparsityPattern.o OpenSMOKELoad.o -o lib/libOpenSMOKEppMath.so 

ar rcs lib/libOpenSMOKEppMath.a OpenSMOKEBaseClass.o OpenSMOKEClass.o OpenSMOKEStdInclude.o OpenSMOKEFunctions.o Conversions.o OpenSMOKEUtilities.o OpenSMOKEVector.o OpenSMOKEMatrix.o OpenSMOKEBandMatrix.o OpenSMOKE_MatrixSparsityPattern.o OpenSMOKELoad.o

rm -f OpenSMOKEBaseClass.o OpenSMOKEClass.o OpenSMOKEStdInclude.o OpenSMOKEFunctions.o Conversions.o OpenSMOKEUtilities.o OpenSMOKEVector.o OpenSMOKEMatrix.o OpenSMOKEBandMatrix.o OpenSMOKE_MatrixSparsityPattern.o OpenSMOKELoad.o


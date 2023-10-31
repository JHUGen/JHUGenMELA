{
  TString LIBCOLLIER = "libcollier.so";
  TString LIBMCFM = "libmcfm_709.so";
  TString LIBJHUGENMELA = "libjhugenmela.so";
  TString LIBMELA = "libJHUGenMELAMELA.so";

  TString LIBPATH = "${MELA_LIB_PATH}/";

  TString LIBMELADIR = LIBPATH;
  if (gSystem->FindDynamicLibrary(LIBMELA)) LIBMELADIR = "";

  TString LIBCOLLIERDIR = LIBPATH;
  if (gSystem->FindDynamicLibrary(LIBCOLLIER)) LIBCOLLIERDIR = "";

  TString LIBMCFMDIR = LIBPATH;
  if (gSystem->FindDynamicLibrary(LIBMCFM)) LIBMCFMDIR = "";

  gInterpreter->AddIncludePath("${ROOFITSYS}/include/");
  gInterpreter->AddIncludePath(LIBPATH+"../../interface/");
  //////////////////////////////////////
  //these explicit loads are required on
  //some machines but not others
  //not entirely sure why
  //either way, they shouldn't hurt
  gSystem->Load("libRooFit");
  gSystem->Load("libPhysics");
  gSystem->Load("libgfortran");
  //////////////////////////////////////
  gSystem->Load(LIBCOLLIERDIR + LIBCOLLIER);
  gSystem->Load(LIBMELADIR + LIBMELA);
  gSystem->Load(LIBMCFMDIR + LIBMCFM);
}

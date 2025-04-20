#include "AtlasStyle.C"
#include "myAtlasUtils.C"
#include "MakeFit_Bins.C"

void rootlogon()
{
  // Load ATLAS style
  //gROOT->LoadMacro("AtlasStyle.C");
  //gROOT->LoadMacro("myAtlasUtils.C");
  SetAtlasStyle();
}
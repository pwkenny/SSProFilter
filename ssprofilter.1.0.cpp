//  SSProFilter Version 1.0
//  Author: Peter W Kenny
//  This source code is provided as supporting information for:
//  PW Kenny, CA Montanari, IM Prokopcyk (2013) ClogPalk: A method for predicting alkane/water
//  partition coefficient. JCAMD 27:*-*  (If citing please consult journal for page numbers 




#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include "openeye.h"
#include "oechem.h"

using namespace OEChem;
using namespace OESystem;
using namespace std;

const char *InterfaceData = 
"!BRIEF UsingOEInterfaceHelp [-o] <output> [-i] <input> [-s] <smarts> [-v] <vbind> \n"
"!PARAMETER -input\n"
"  !ALIAS -i\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF Input file\n"
"  !DETAIL\n"
"Input file of molecules\n"
"!END\n"
"!PARAMETER -output\n"
"  !ALIAS -o\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF The count file\n"
"  !DETAIL\n"
"Output file \n"
"!END\n"
"!PARAMETER -smarts\n"
"  !ALIAS -s\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF SMARTS\n"
"  !DETAIL\n"
"SMARTS definition file\n"
"Lines in this file starting with # will be treated as comments\n"
"Different formats of SMARTS file are used for profiling and filtering\n"
"Format for profiling: SMARTS.string   SMARTS.string\n"
"Format for filtering: SMARTS.string   minimum.matches   maximum.matches\n"
"!END\n"
"!PARAMETER -vectorbind\n"
"  !ALIAS -v\n"
"  !TYPE string\n"
"  !REQUIRED false\n"
"  !BRIEF Vector bindings\n"
"  !DETAIL\n"
"Vector binding definition file (optional)\n"
"Lines in this file starting with # will be treated as comments\n"
"Format: vector.binding.name vector.binding.definition\n"
"!END\n"
"!PARAMETER -type\n"
"  !ALIAS -t\n"
"  !TYPE string\n"
"  !REQUIRED true\n"
"  !BRIEF Calculation type\n"
"  !DETAIL\n"
"Sets calculation type (filter or profile) \n"
"!END\n" ;



int main(int argc , char** argv)

{

 OEInterface itf(InterfaceData, argc, argv);

 OEGraphMol mol;
 OESubSearch *ss;

 string smiles , strname;
 string *smt_string , *smt_name;
 string vbname , vbdef , line , calctype ;
 vector< pair<string, string> > definitions;
 unsigned int *i_low;
 unsigned int *i_high;
 int i , i_smt , i_accept;
 unsigned int count;

 oemolistream infil;
 oemolostream outfil;
 ifstream smtfil;
 ifstream vbfil;
 ofstream countfil;



 calctype = itf.Get<std::string>("-t");
 infil.open(itf.Get<std::string>("-i"));
 smtfil.open(itf.Get<std::string>("-s").c_str());

  if (itf.Has<string>("-v"))
   {
     vbfil.open(itf.Get<std::string>("-v").c_str());
     if (!vbfil) OEThrow.Error("Unable to open vector binding file\n",itf.Get<std::string>("-v").c_str());
   }

  if ( calctype == "profile" )
   {
     cout << "check count\n";
     countfil.open(itf.Get<std::string>("-o").c_str());
     if (!countfil) OEThrow.Error("Unable to open profile file\n",itf.Get<std::string>("-o").c_str());
   }
  else 
   {
     outfil.open(itf.Get<std::string>("-o"));
     if (!outfil) OEThrow.Error("Unable to open output file\n",itf.Get<std::string>("-o").c_str());
   }

  if (!infil) OEThrow.Error("Unable to open input file\n",itf.Get<std::string>("-i").c_str());
  if (!smtfil) OEThrow.Error("Unable to open SMARTS file\n",itf.Get<std::string>("-s").c_str());

  // Read vector bindings

  if (vbfil)
    {
      while(getline(vbfil,line))
        {
          if (!vbfil.eof())
            {  
               if ( line.find("#") != 0) 
                 {
                    istringstream iss(line);
                    iss >> vbname >> vbdef; 
                    definitions.push_back(pair<string,string>( vbname, vbdef ));
                 }
            }
        }
    }
//  Count substructural definitions
 
 i_smt = 0;
 while(getline(smtfil,line))
   {
     if (!smtfil.eof())
        {  
           if ( line.find("#") != 0) 
              {
                  i_smt++;
               }
         }
   }

 // Prepare the arrays and read substructural definitions

 smtfil.close();
 smtfil.open(itf.Get<std::string>("-s").c_str());
 ss = new OESubSearch [i_smt];
 smt_name = new string [i_smt];
 smt_string = new string [i_smt];
 i_low = new unsigned int [i_smt];
 i_high = new unsigned int [i_smt]; 
 i_smt = 0;
 while(getline(smtfil,line))
   {
     if (!smtfil.eof())
        {  
           if ( line.find("#") != 0) 
              {
                  cout << "check line" << line << "\n";
                  istringstream iss(line);
                  if ( calctype == "filter" )                  
                     iss >> smt_string[i_smt] >> i_low[i_smt] >> i_high[i_smt];
                  if ( calctype == "profile" )
		    iss >> smt_name[i_smt] >> smt_string[i_smt];
                  OESmartsLexReplace(smt_string[i_smt], definitions);
                  if (!ss[i_smt].Init(smt_string[i_smt].c_str()))
                      OEThrow.Error("Unable to parse %s",smt_string[i_smt].c_str());
                  i_smt++;
               }
         }
   }

 if ( calctype == "profile" )
   {
     countfil << "Name" ;
     for ( i= 0 ; i < i_smt ; i++ ) 
         countfil << " " << smt_name[i];
     countfil << "\n";
   }
   


// Read input molecules

 while (OEReadMolecule(infil, mol))
   {
     if (!infil.eof()) {
        i_accept = 1;
        i = 0;
        if ( calctype == "profile" )  countfil << mol.GetTitle();
        while( ( i < i_smt ) && ( i_accept == 1) )
	//        for ( i = 0 ; i < i_smt ; i++) 
            {
              count = 0;
              for (OEIter<OEMatchBase> match = ss[i].Match(mol); match; ++match, ++count)
                 { }
              if ( calctype == "profile" )
		{ 
                  countfil << " " << count;
		 } 
              else if ( count < i_low[i] || count > i_high[i] )
	         {
                   i_accept = 0;
                 }
              i++;
	    }
        if ( calctype == "profile" ) 
	  {
            countfil << "\n";
          }
        else
          {
            if ( i_accept == 1) OEWriteMolecule(outfil, mol);
          }
        mol.Clear();
   }
 }
 infil.close(); 

}








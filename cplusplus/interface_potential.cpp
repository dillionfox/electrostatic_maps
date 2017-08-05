#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <ctime>
#include <malloc.h>
#include </home/dillion/pkg/surf-master/depends/xdrfile-1.1.4/include/xdrfile_xtc.h> // xdr include file

using namespace std;
//#include "calphi.h"
//#include "array_3D.h"
//#include "readpro.h"
//#include "Marching_Cube.h"

//***********************************************************************
// ::Version 2::
// Dillion Fox, 2017
//
// Original code written by R. Remsing, 2014
//
// This code will calculate the long-ranged component of the
// electrostatic potential arising from *WATER ONLY* at each grid point
// along the Willard-Chandler instantaneous interface
// (Willard, Chandler, J. Phys. Chem. B 2010, 114, 1954).
//
// For details of the method, see:
// Remsing, Weeks, J. Phys. Chem. B 2015, 119, 9268-9277
//
// The code will also compute the distribution of the potential
// over all grid points.
//***********************************************************************

// To compile:
// ./compile_interface_potential.sh
//
// To run:
// ./interface_potential.out protII.xtc trajectory.xtc
//

// As mentioned in the instantaneous interface (II) code comments, we assume the II is for a rigid
// protein, and so only one interface configuration is needed. The code will need to be modified
// to incorporate many configurations of the interface if desired.

int main(int argc, char **argv) {
  /* XTC variables */
  XDRFILE *xdin;                // the xdr file holder
  XDRFILE *xdout;
  XDRFILE *intin;               // xdr file holder for interface coordinates

  // USER INPUT
  int npro = 14;                // number of protein atoms
  int lastfr = 6;               // last frame of the trajectory to analyze
  double nconf = 6.0;           // total number of configurations in the trajectory
  int first_water = 14;         // first water index
  int last_water = 16645;       // last water index

  int natoms=0;                 // number of atoms
  int nintatoms=0;              // number of interface atoms
  int step;                     // the step counter for the simulation
  float time;                   // simulation time
  matrix box;                   // box coordinates in a 3x3 matrix
  matrix intbox;                // box coords for interface file
  rvec *xtc_coor, *int_coor;    // atom coordinates in a 2-D matrix
  float prec;                   // precision of the xtc file
  float intprec;

  /* Program variables */
  int frame = 0;
  int cutfr = 0;                // skip <cutfr> frames at the beginning of the trajectory before analyzing
  int countfr = 0;
  int skipfr = 1;               // skip every <skipfr> frames


  float xij,yij,zij,rij;
  float xiji,yiji,ziji,riji;
  
  float eps0;
  float sigma=0.45;             // Smoothing length for Local Molecualr Field (LMF) theory decomposition of the electrostatic potential
  float chg[3];                 // Charges of the water oxygen and hydrogen atoms *(code is specific for 3-site water models)*
                                // Current charges are for SPC/E water
  chg[0] = -0.834;
  chg[1] = 0.417;
  chg[2] = chg[1];
  
  /* Memory allocation to read coordinates */
  // This is for the interface
  read_xtc_natoms(argv[1], &nintatoms);
  int_coor = (rvec *) malloc(nintatoms*sizeof(rvec));
  if(int_coor==0){
    cout<<"Insufficient memory to load .xtc file.\n";
    return 0;
  }

  float VRS[nintatoms];
  float VRSt[nintatoms];
  for(int i=0;i<nintatoms;i++)
  {
    VRS[i] = 0.0;
    VRSt[i] = 0.0;
  }

  /* Memory allocation to read coordinates - This is for the system */
  read_xtc_natoms(argv[2], &natoms);
  xtc_coor = (rvec *) malloc(natoms*sizeof(rvec));
  if(xtc_coor==0)
  {
    cout << "Insufficient memory to load .xtc file.\n";
    return 0;
  }

  /* Open the xtc file and loop through each frame. */
  xdin=xdrfile_open(argv[2],"r");
  
  intin = xdrfile_open(argv[1],"r"); // interface coords

  read_xtc(intin, nintatoms, &step,&time,intbox,int_coor,&intprec);
  while( ! read_xtc(xdin, natoms, &step, &time, box, xtc_coor, &prec) ) 
  {
    frame ++;   
    if (frame > cutfr && frame <= lastfr && (frame-1)%skipfr == 0 )
    {
      cout << "Analyzing frame " << frame << endl;
      countfr ++;
      for(int ii=0; ii<nintatoms; ii++) // Loop over interface points
      {
        // Compute the LR potential at the interface site due to the waters
        //for( int i=npro; i<natoms; i+=3 )
        for( int i=first_water; i<last_water; i+=3 )
        {
          for( int j=0;j<3;j++)
          {
      	    xij = xtc_coor[i+j][0] - int_coor[ii][0]; // water_coor - interface_coor
      	    yij = xtc_coor[i+j][1] - int_coor[ii][1];
      	    zij = xtc_coor[i+j][2] - int_coor[ii][2];
      	    rij = sqrt( xij*xij + yij*yij + zij*zij );
      	    
      	    xiji = xij - int( xij/box[0][0] + 0.5 )*box[0][0];
      	    yiji = yij - int( yij/box[1][1] + 0.5 )*box[1][1];
      	    ziji = zij - int( zij/box[2][2] + 0.5 )*box[2][2];
      	    riji = sqrt( xiji*xiji + yiji*yiji + ziji*ziji );
      	    
      	    // Still need to add the 4*pi*eps0 in the denominator if you want SI units!!!!
      	    VRS[ii] += chg[j] * erf(rij/sigma) / (rij);
          }
        }
        //cout << "natoms: " << natoms << ", VRS[ii]: " << VRS[ii] << ", erf(): " << erf(rij/sigma) << ", rij: " << rij << endl;
      }
      if( countfr%1 == 0 )
      {
        for( int ii=0; ii<nintatoms; ii++)     // Now we have the potential at each grid point.
          VRSt[ii] = VRS[ii]/float(countfr);   // Let's compute the distribution of the potential
                                               // at all grid points except (0,0,0)
        float minVRS = 10000000.0;
        float maxVRS = -10000000.0;
        
        // First find min and max of VRS
        for( int ii=0; ii<nintatoms; ii++)
        {
          if(int_coor[ii][0] > 0.0)
          {
            if(VRSt[ii] < minVRS)
              minVRS = VRSt[ii];
            if(VRSt[ii] > maxVRS)
              maxVRS = VRSt[ii];
          }
        }
      
        int nbinV = 50; // Set up bins
        float dV = (maxVRS-minVRS)/float(nbinV);
        float PV[nbinV];
        float refVal;
        int iV;
        cout << "Calculating Distribution..." << endl;
        
        for( int ii=0;ii<nbinV;ii++)
        	PV[ii] = 0.0;
        
        for(int ii=0;ii<nintatoms;ii++)
        {
          //cout << "int_coor[ii][0]: " << int_coor[ii][0] << endl;
          if(int_coor[ii][0] > 0.0)
          {
            iV = int( (VRSt[ii] - minVRS)/dV );
            PV[iV] += 1.0;
          }
          else
            refVal = VRSt[ii];
        }
        cout << "REFVAL: " << refVal << endl;
        
        float sum=0.0;
        for(int ii=0;ii<nbinV;ii++)
          sum += PV[ii];
        
        ofstream pfile;
        pfile.open("VDist.dat");
        float VVV;
        for(int ii=0;ii<nbinV;ii++)
        {
          PV[ii]/=sum;
          VVV = minVRS+ii*dV;
          pfile << "\t" << VVV << "\t" << PV[ii] << "\n";
        }
        pfile.close();
        
        ofstream vtfFile; // Write out a VTF file for viewing in VMD
        ofstream vsfFile;
        ofstream vtfS;

        string filename;
        string head = "Map";
        string ext = ".vtf";
        stringstream ss;
        ss << time;
        filename = head + ss.str() + ext;

        vtfFile.open(filename.c_str());
        vtfS.open("Map-scaled.vtf");
        vsfFile.open("Map.vsf");
        vtfFile << "# Structure Block" << endl;
        vsfFile << "# Structure file" << endl;
        vtfFile << "atom default radius 1.0 name N" << endl;
        vsfFile << "atom default radius 1.0 name N" << endl;
        vtfS << "# Structure Block" << endl;
        vtfS << "atom default radius 1.0 name N" << endl;

        for(int ii=0;ii<nintatoms;ii++)
        {
          if( int_coor[ii][0] > 0.0)
          {
            vtfFile << " atom " << ii << " b " << VRSt[ii]-refVal << endl;
            vsfFile << " atom " << ii << " b " << VRSt[ii]-refVal << endl;
            vtfS << " atom " << ii << " b " << ( VRSt[ii]-refVal - minVRS + refVal ) / ( maxVRS-refVal - minVRS + refVal ) << endl;
          }
        }
        vtfFile << "# timestep block" << endl;
        vtfFile << "timestep indexed" << endl;
        vtfFile << " pbc  " << intbox[0][0]*10.0 << " " << intbox[1][1]*10.0 << " " << intbox[2][2]*10.0 << endl;
        vtfS << "# timestep block" << endl;
        vtfS << "timestep indexed" << endl;
        vtfS << " pbc  " << intbox[0][0]*10.0 << " " << intbox[1][1]*10.0 << " " << intbox[2][2]*10.0 << endl;
        
        for(int ii=0;ii<nintatoms;ii++)
        {
          if( int_coor[ii][0] > 0.0 )
          {
            vtfFile << " " << ii << " " << int_coor[ii][0]*10.0 << " " << int_coor[ii][1]*10.0 << " " << int_coor[ii][2]*10.0 << endl;
            vtfS << " " << ii << " " << int_coor[ii][0]*10.0 << " " << int_coor[ii][1]*10.0 << " " << int_coor[ii][2]*10.0 << endl;
          }
        }
        vtfFile.close();
        vsfFile.close();
      } // End if statement for periodic printing the output
    }
  }
  
  for( int ii=0; ii<nintatoms; ii++) // Print final averages
    VRS[ii] /= float(countfr);

  float minVRS = 10000000.0;  // Now we have the potential at each grid point.
  float maxVRS = -10000000.0; // Let's compute the distribution of the potential
                              // at all grid points except (0,0,0)

  for( int ii=0; ii<nintatoms; ii++) // First find min and max of VRS
  {
    if(int_coor[ii][0] > 0.0)
    {
      if(VRS[ii] < minVRS)
          minVRS = VRS[ii];
      if(VRS[ii] > maxVRS)
          maxVRS = VRS[ii];
    }
  }
  
  int nbinV = 50; // Set up bins
  float dV = (maxVRS-minVRS)/float(nbinV);
  float PV[nbinV];
  float refVal;
  int iV;
    
  // Calculate distribution of VR
  cout << "Calculating Distribution..." << endl;

  for( int ii=0;ii<nbinV;ii++)
    PV[ii] = 0.0;

  for(int ii=0;ii<nintatoms;ii++)
  {
    if(int_coor[ii][0] > 0.0)
    {
      iV = int( (VRS[ii] - minVRS)/dV );
      PV[iV] += 1.0;
    }
    else
        refVal = VRS[ii];
  }
  cout << "DENOMINATOR: " << (maxVRS -refVal - minVRS + refVal) << ", REFVAL: " << refVal <<  endl;

  float sum=0.0;
  for(int ii=0;ii<nbinV;ii++)
    sum += PV[ii];

  ofstream pfile; // Write out distribution of VR
  pfile.open("VDist.dat");
  float VVV;
  for(int ii=0;ii<nbinV;ii++)
  {
    PV[ii]/=sum;
    VVV = minVRS+ii*dV;
    pfile << "\t" << VVV << "\t" << PV[ii] << "\n";
  }
  pfile.close();

  ofstream vtfFile; // Write out a VTF file for viewing in VMD
  ofstream vsfFile;
  ofstream vtfS;
  vtfFile.open("Map.vtf");          // VTF file of the map
  vtfS.open("Map-scaled.vtf");      // VTF file of the map scaled
  vsfFile.open("Map.vsf");          // VSF file of the map

  vtfFile << "# Structure Block" << endl;
  vsfFile << "# Structure file" << endl;
  vtfFile << "atom default radius 1.0 name N" << endl;
  vsfFile << "atom default radius 1.0 name N" << endl;

  vtfS << "# Structure Block" << endl;
  vtfS << "atom default radius 1.0 name N" << endl;
  for(int ii=0;ii<nintatoms;ii++)
  {
    if( int_coor[ii][0] > 0.0)
    {
      vtfFile << " atom " << ii << " b " << VRS[ii]-refVal << endl;
      vsfFile << " atom " << ii << " b " << VRS[ii]-refVal << endl;
      vtfS << " atom " << ii << " b " << ( VRS[ii]-refVal - minVRS + refVal ) / ( maxVRS-refVal - minVRS + refVal ) << endl;
    }
  }
  vtfFile << "# timestep block" << endl;
  vtfFile << "timestep indexed" << endl;
  vtfFile << " pbc  " << intbox[0][0]*10.0 << " " << intbox[1][1]*10.0 << " " << intbox[2][2]*10.0 << endl;
  
  vtfS << "# timestep block" << endl;
  vtfS << "timestep indexed" << endl;
  vtfS << " pbc  " << intbox[0][0]*10.0 << " " << intbox[1][1]*10.0 << " " << intbox[2][2]*10.0 << endl;

  for(int ii=0;ii<nintatoms;ii++)
  {
    if( int_coor[ii][0] > 0.0 )
    {
      vtfFile << " " << ii << " " << int_coor[ii][0]*10.0 << " " << int_coor[ii][1]*10.0 << " " << int_coor[ii][2]*10.0 << endl;
      vtfS << " " << ii << " " << int_coor[ii][0]*10.0 << " " << int_coor[ii][1]*10.0 << " " << int_coor[ii][2]*10.0 << endl;
    }
  }
  vtfFile.close();
  vsfFile.close();

  xdrfile_close(xdin);
  xdrfile_close(intin);

  return 0;
}

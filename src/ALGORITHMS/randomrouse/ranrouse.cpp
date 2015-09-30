#include <iostream>
#include <fstream>
#include <string>
using namespace std;

bool root;

#include "ap.cpp"

#include "blas.cpp"
#include "rotations.cpp"
#include "tdevd.cpp"
#include "sblas.cpp"
#include "reflections.cpp"
#include "tridiagonal.cpp"
#include "sevd.cpp"

#include "readdata.h"

#include "mpi.h"

#include "povray.h"

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

void print_matrix(ap::real_2d_array& m, int n) {
  int i,j;
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      cout << " " << m(i,j); // first index:row
    }
    cout << endl;
  }
}

void print_vector(ap::real_1d_array& v, int n) {
  int i;
  for(i=0;i<n;i++) {
    cout << " " << v(i) << endl;
  }
}

void sort_vector(ap::real_1d_array& v, int n) {
  int i, j=n-1;
  double temp;
  bool changed;
  do {
    changed=false;
    for(i=0;i<j;i++) {
      if(v(i)>v(i+1)) {
        SWAP(v(i), v(i+1));
        changed=true;
      }
    }
    j--;
  } while(changed && (j>0));
}

double get_random_lognorm(double tau_m, double dtau, double lower, double upper) {
  double x,y;
  do
  {
    x = lower + (upper-lower)*ap::randomreal();
    y = ap::randomreal();
  } while (y > exp(-((log(x)-log(tau_m))/dtau)*((log(x)-log(tau_m))/dtau)));
  return x;
}  

int main(int argc, char *argv[]) {
  if(root) {
    cout << "Current time: " << std::flush;
    system("date");
    cout << endl;
    cout << "Brodecks Random Rouse\n";
    cout << "This is v0.6 of February 26th 2010\n";
    cout << "\n";
    cout << "This program is the result of INTENSIVE testing\n";
    cout << "and comparison with very simple systems and\n";
    cout << "assumptions! Please report any bugs you encounter\n";
    cout << "to brodeck@gmail.com\n";
    cout << "\n";
    cout << "Hints:\n";
    cout << "How to get zeta_0 from tau_R in the pure Rouse?\n";
    cout << "-> tau_R=N^2/pi^2*zeta_0\n";
    cout << "\n";
    cout << "This way you can connect this program to real exp. values\n";
    cout << "\n";
    cout << "-------------------------------------------------------\n";
    cout << "\n";
  }

  if(argc<=12) { // n-1 parameters needed
    cout << "Need parameters:" << endl;
    //cout << "  ranrouse.exe chains MSD_entries N(monomers) zeta_0 l t_min t_max errorcheck(1/0) debug(1/0) povray(1/0) prefix friction_dist(use -1 for help) additional_parameters(see help)" << endl;
    //cout << "  ranrouse.exe 10 150 43 1.0 7.24 0.1 1000 1 0 run2 0" << endl;
    cout << "  ranrouse.exe chains MSD_entries N(monomers) Wl4 Ree t_min t_max errorcheck(1/0) debug(1/0) povray(1/0) prefix friction_dist(use -1 for help) additional_parameters(see help)" << endl;
    cout << "  ranrouse.exe 10 150 43 15091 135 0.1 1000 1 0 run2_ 0" << endl;
    cout << "tau_min and tau_max multiples of the Rouse time" << endl;
    cout << "\n";
    cout << "Parameter description:\n";
    cout << "chains:      the number of chains that are simulated (important for random parameters)\n";
    cout << "             if there is a fixed friction coefficients one chain is sufficient\n";
    cout << "MSD_entries: number of MSD (and S_coherent) entries in the final functions (t-values)\n";
    cout << "N:           number of beads in the chain\n";
    cout << "Wl4:         Rouse rate W times l^4 (bond length), determined experimentally\n";
    cout << "Ree:         end-to-end distance of the chain (this defines l)\n";
    cout << "t_min:       smallest time calculated\n";
    cout << "t_max:       longest time calculated\n";
    cout << "errorcheck:  if set to one will check the consistency of the code using different methods\n";
    cout << "             of calculation. Takes much longer!\n";
    cout << "debug:       display (many many) debug messages\n";
    cout << "povray:      output povray render files\n";
    cout << "prefix:      prefix for the file output\n";
    cout << "friction_..: type of friction distribution. set to -1 for options!\n";
    cout << "add_param:   additional parameters for the friction distribution\n";
    cout << "\n";
    return -1;
  }

  // MPI BEGIN
  int pid, messagetag=99, numprocs, namelen;
  int root_i=0;

  root=false;

  /// This code can run in parallel, so lets setup MPI
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processor_name, &namelen);
  cout << " --> Dienstantrittsmeldung: Hi! I am processor " << pid << " on " << processor_name << " out of " << numprocs << endl << flush;

  if(pid==root_i) root=true;

  srand((unsigned)time(NULL)+pid+rand());
  // MPI END

  int i,j,x,y,p,q;

  cout.precision(5);
  cout.setf(ios::fixed,ios::floatfield);

  // get parameters from command line input
  
  int chains=(atoi(argv[1]));
  int MSD_entries=(atoi(argv[2]));
  int N=(atoi(argv[3]));
  //double zeta_0=((double)atof(argv[4]));
  //double l=((double)atof(argv[5]));

  double Wl4=((double)atof(argv[4]));
  double Ree=((double)atof(argv[5]));

  double l=Ree/sqrt((double)N);
  double zeta_0=Ree*Ree*Ree*Ree/Wl4/(double)(N*N);

  double tau_R=zeta_0*(double)(N*N)/(M_PI*M_PI);

  double t_min=((double)atof(argv[6]))*tau_R/(N*N); // 0.01*
  double t_max=((double)atof(argv[7]))*tau_R; // 1000*

  bool debug=false;
  if(atoi(argv[9])==1) debug=true;

  bool echeck=debug;
  if(atoi(argv[8])==1) echeck=true;

  bool povray=false;
  if(atoi(argv[10])==1) povray=true;

  string prefix=argv[11];

  int friction_dist=(atoi(argv[12]));
  if(friction_dist==-1) {
    if(root) {
      cout << "Possible paramters for friction distributions\n";
      cout << "\n";
      cout << "-1 - Display this help\n";
      cout << " 0 - pure Rouse, constant frictions\n";
      cout << " 1 - Fix one bead\n";
      cout << "     Additional parameter 1: Friction factor of fixed bead (e.g. 10e10)\n";
      cout << "     Additional parameter 2: index of fixed bead (starting at 0)\n";
      cout << " 2 - Fix several beads at random position\n";
      cout << "     Additional parameter 1: Friction factor of fixed bead (e.g. 10e10)\n";
      cout << "     Additional parameter 2: number of fixed beads\n";
      cout << " 3 - Fix one bead, average over whole chain\n";
      cout << "     Additional parameter 1: Friction factor of fixed bead (e.g. 10e10)\n";
      cout << "     Parameter \"number of chains\" will be overwritten!\n";
      cout << " 4 - Read friction coefficients from file\n";
      cout << "     Additional parameter 1: filename\n";
      cout << "     Parameters \"number of chains\", \"Wl4\", will be overwritten!\n";
      cout << "     Data format: first column index, second column friction.\n";
      cout << "     Always reads N lines, then starts next chain.\n";
      cout << "\n";
      cout << "\n";
    }
    return 0;
  }

  double dextra_par1;
  int iextra_par1;
  string zeta_dist_filename;

  if(friction_dist==1)  {
    if(argc<=14) {
      cout << "Not enough parameters for this distribution!\n";
      cout << "Exiting... see help (-1) for details!\n";
      return -1;
    }
    dextra_par1=((double)atof(argv[13]));
    iextra_par1=(atoi(argv[14]));
  }
  if(friction_dist==2)  {
    if(argc<=14) {
      cout << "Not enough parameters for this distribution!\n";
      cout << "Exiting... see help (-1) for details!\n";
      return -1;
    }
    dextra_par1=((double)atof(argv[13]));
    iextra_par1=(atoi(argv[14]));
    if(iextra_par1>=N) {
      cout << "Error, can not have more fixed beads than N-1" << endl;
      return -1;
    }
  }
  if(friction_dist==3)  {
    if(argc<=13) {
      cout << "Not enough parameters for this distribution!\n";
      cout << "Exiting... see help (-1) for details!\n";
      return -1;
    }
    dextra_par1=((double)atof(argv[13]));
    chains=N;
  }
  double *inputfrict=NULL;
  int total_cols;
  if(friction_dist==4)  {
    if(argc<=12) {
      cout << "Not enough parameters for this distribution!\n";
      cout << "Exiting... see help (-1) for details!\n";
      return -1;
    }
    zeta_dist_filename=argv[13];

    i=read_column(zeta_dist_filename, inputfrict, 1, 1); // first column, ignore first line
    cout << i << " lines found in file!" << endl;
    if(i%N!=0) cout << "ERROR: Does not compute with your value of N: " << N << endl;
    chains=i/N;
    cout << chains << " chains!" << endl;

    Ree=0.0;
    l=0.0;
    tau_R=0.0;
    zeta_0=0.0;
  }

  int chain;

  //prefix+="_";

  double ERROR_MARGIN=1e-6;

  //double tau_R=tau*(double)(N*N)/(M_PI*M_PI);

  double D;

  ap::real_1d_array msd;
  msd.setbounds(0, MSD_entries-1);
  ap::real_1d_array msd_woD;
  msd_woD.setbounds(0, MSD_entries-1);

  ap::real_1d_array times;
  times.setbounds(0, MSD_entries-1);

  int numQ=41;

  // list of Q-values that will be calculated!
  double Qval[numQ];
  Qval[0]=0.045;
  Qval[1]=0.050;
  Qval[2]=0.055;
  Qval[3]=0.060;
  Qval[4]=0.065;
  Qval[5]=0.070;
  Qval[6]=0.075;
  Qval[7]=0.080;
  Qval[8]=0.085;
  Qval[9]=0.090;
  Qval[10]=0.095;
  Qval[11]=0.100;
  Qval[12]=0.105;
  Qval[13]=0.110;
  Qval[14]=0.115;
  Qval[15]=0.120;
  Qval[16]=0.125;
  Qval[17]=0.130;
  Qval[18]=0.135;
  Qval[19]=0.135;
  Qval[20]=0.140;
  Qval[21]=0.145;
  Qval[22]=0.150;
  Qval[23]=0.175;
  Qval[24]=0.20;
  Qval[25]=0.25;
  Qval[26]=0.30;
  Qval[27]=0.40;
  Qval[28]=0.50;
  Qval[29]=0.60;
  Qval[30]=0.70;
  Qval[31]=0.80;
  Qval[32]=0.90;
  Qval[33]=1.00;
  Qval[34]=1.10;
  Qval[35]=1.20;
  Qval[36]=1.30;
  Qval[37]=1.40;
  Qval[38]=1.50;
  Qval[39]=1.60;
  Qval[40]=1.70;

  ap::real_2d_array Scoh; // coherent scattering functions
  Scoh.setbounds(0, numQ, 0, MSD_entries);

  ap::real_2d_array Sinc; // coherent scattering functions
  Sinc.setbounds(0, numQ, 0, MSD_entries);

  ap::real_2d_array Scoh_woD; // coherent scattering functions without diffusion
  Scoh_woD.setbounds(0, numQ, 0, MSD_entries);

  for(i=0;i<numQ;i++) {
    for(j=0;j<MSD_entries;j++) {
      Sinc(i, j)=0.0;

      Scoh(i, j)=0.0;
      Scoh_woD(i, j)=0.0;
    }
  }

  double dt=(log10(t_max)-log10(t_min))/(double)MSD_entries;

  for(i=0;i<MSD_entries;i++) {
    msd(i)=0.0;
    msd_woD(i)=0.0;
    times(i)=pow(10, log10(t_min)+i*dt);
  }

  double sums;

  MPI_Barrier(MPI_COMM_WORLD);

  if(root) {
    cout << "Current time: " << std::flush;
    system("date");
    cout << endl;
    cout << "Brodecks Random Rouse\n";
    cout << "This is v0.6 of February 26th 2010\n";
    cout << "\n";
    cout << "This program is the result of INTENSIVE testing\n";
    cout << "and comparison with very simple systems and\n";
    cout << "assumptions! Please report any bugs that you MIGHT\n";
    cout << "encounter (even though I am pretty sure there are none)\n";
    cout << "to brodeck@gmail.com\n";
    cout << "\n";
    cout << "Hints:\n";
    cout << "How to get zeta_0 from tau_R?\n";
    cout << "Easy: tau_R=N^2/pi^2*zeta_0\n";
    cout << "\n";
    cout << "This way you can connect this program to real exp. values\n";
    cout << "\n";
    cout << "-------------------------------------------------------\n";
    cout << "\n";
    cout << "Random number test, seed implemented: " << ap::randomreal() << endl;
    cout << "Parameters: " << endl;

    cout << "N: " << N << endl;
    cout << "chains: " << chains << endl;
    cout << "l: " << l << endl;
    cout << "zeta_0: " << zeta_0 << endl;
    cout << "tau_R: " << tau_R << endl;
    cout << "Wl^4: " << l*l*l*l/(zeta_0) << endl;
    cout << "msd t_max set to: " << t_max << endl;
    cout << "msd t_min set to: " << t_min << endl;

    cout << "Error check: ";
    if(echeck) cout << "YES"; else cout << "NO";
    cout << endl;
    cout << "Debug: ";
    if(debug) cout << "YES"; else cout << "NO";
    cout << endl;
  }

  ap::real_1d_array friction;
  friction.setbounds(0, N-1);

  double friction_av;

  ap::real_1d_array mobility;
  mobility.setbounds(0, N-1);

  ap::real_1d_array A_diag; // mobility matrix diagonal
  A_diag.setbounds(0, N-2);
  ap::real_1d_array A_tridiag; // mobility matrix tridiagonal
  A_tridiag.setbounds(0, N-3);

  ap::real_1d_array er_av; // eigenvalues averaged
  er_av.setbounds(0, N-2);
  for(i=0;i<N-1;i++) {
    er_av(i)=0;
  }

  ap::real_2d_array ev; // eigenvectors
  ev.setbounds(0, N-2, 0, N-2);

  ap::real_2d_array c; // needed for msd, c(i,j), see document
  c.setbounds(0, N-1, 0, N-2);

  ap::real_2d_array phi; // phi_ij correlators
  phi.setbounds(0, N-1, 0, N-1);

  for(chain=0;chain<chains;chain++) {
    if(chain%numprocs!=pid) continue;

    cout << "pid " << pid << " - Working on chain " << chain << " of " << chains << endl;

    friction_av=0;

    if(debug) cout << "Frictions: " << endl;
    for(i=0;i<N;i++) {
      // Set friction coefficients
      friction(i)=zeta_0;
    }

    if(friction_dist==1)  { // fix one bead
      friction(iextra_par1)*=dextra_par1;
    }
    if(friction_dist==2)  { // fix several beads
      for(i=0;i<iextra_par1;i++) {
        do {
          j=ap::randominteger(N);
          if(friction(j)-zeta_0<ERROR_MARGIN) {
            friction(j)*=dextra_par1;
            break;
          }
        } while(true);
      }
    }
    if(friction_dist==3)  { // fix bead with index chain
      friction(chain)*=dextra_par1;
    }
    if(friction_dist==4)  { // read frictions from file
      for(i=0;i<N;i++) {
        friction(i)=1/inputfrict[chain*N+i];
        //friction(i)
        //friction(i)=64.10705725;
      }
    }

    for(i=0;i<N;i++) {
      mobility(i)=1/friction(i);
      friction_av+=friction(i);
      if(debug) cout << "  friction " << i << " : " << friction(i) << "\t mobility : " << mobility(i) << endl;
    }

    D=l*l/(3.0*friction_av);
    friction_av/=double(N);

    if(debug) cout << "Diffusion: " << D << endl;

    for(i=0;i<N-1;i++) {
      A_diag(i)=(mobility(i)+mobility(i+1));
    }
    for(i=0;i<N-2;i++) {
      A_tridiag(i)=-mobility(i+1);
    }

    if(debug) {
      cout << "Diagonal and tridiagonal vectors: " << endl;
      print_vector(A_diag, N-1);
      cout << endl;
      print_vector(A_tridiag, N-2);
      cout << endl;
    }

    ap::real_2d_array tmpm; // eigenvectors
    if(debug) {
      tmpm.setbounds(0, N-1, 0, N-1);
      for(i=0;i<N-1;i++) {
        for(j=0;j<N-1;j++) {
          tmpm(i,j)=0.0;
        }
      }
      for(i=0;i<N-1;i++) {
        tmpm(i,i)=A_diag(i);
        if(i<N-2) tmpm(i+1,i)=A_tridiag(i);
        if(i<N-2) tmpm(i,i+1)=A_tridiag(i);
      }
    }

    if(debug) {
      cout << "Mobility matrix:" << endl;
      cout << endl;
      print_matrix(tmpm, N-1);
      cout << endl;
    }

    /*
     * THIS IS THE ALTERNATIVE CALCULATION OF THE EIGENVECTORS
     * However, it yields the same results, of course!
     * It uses a different algorithm...
		if(!smatrixevd(tmpm, N-1, true, 1, A_diag, ev)) { // results now in A_diag
			cout << "Error calculating eigenvectors" << endl;
			return -1;
		}

		cout << endl;
		cout << "EV:" << endl;
		print_vector(A_diag, N-1);
		cout << endl;
		print_matrix(ev, N-1);

		for(i=0;i<N-1;i++) {
			A_diag(i)=(mobility(i)+mobility(i+1));
		}
     */

    // Solve the eigenvalue/eigenvector problem
    if(!smatrixtdevd(A_diag, A_tridiag, N-1, 2, ev)) { // results now in A_diag
      cout << "Error calculating eigenvectors" << endl;
      return -1;
    }

    // Check if eigenvalues are sorted!
    for(i=0;i<N-2;i++) {
      if(A_diag(i)>A_diag(i+1)) {
        cout << "Error sorting your results";
        return -1;
      }
    }

    /* This is no longer necessary! The signs of the eigenvectors do NOT
     * make a difference, as it should be!
		if(debug) cout << "Correcting sign of eigenvectors..." << endl;

		if(debug) {
			for(x=0;x<N-1;x++) {
				for(y=0;y<N-1;y++) {
					sums=0;
					for(i=0;i<N-1;i++) {
						for(j=0;j<N-1;j++) {
							sums+=ev(i,x)*ev(j,y);
						}
					}
					cout << " " << sums;
				}
				cout << endl;
			}
		}
     */

	// output debug information
    if(debug) {
      for(x=0;x<N-1;x++) {
        for(y=0;y<N-1;y++) {
          sums=0;
          for(i=0;i<N-1;i++) {
            for(j=0;j<N-1;j++) {
              sums+=ev(i,x)*ev(j,y);
            }
          }
          if(fabs(sums)<ERROR_MARGIN) cout << " 0";
          else if(sums>0) cout << " +";
          else cout << " -";
        }
        cout << endl;
      }
    }

	// output debug information
    if(debug) {
      cout << "Eigenvalues:" << endl;
      cout << endl;
      print_vector(A_diag, N-1);
      cout << endl;
      cout << "Eigenvectors:" << endl;
      cout << "The first COLUMN corresponds to the first eigenvalue!" << endl;
      cout << endl;
      for(i=0;i<N-1;i++) {
        cout << " " << A_diag(i);
      }
      cout << endl;
      for(i=0;i<N-1;i++) {
        cout << "|" << "-------";
      }
      cout << "|" << endl;
      print_matrix(ev, N-1);
      cout << endl;
    }

    if(debug)
    {
      cout << "Eigenvalues Eigenvector Check" << endl;
    }

	// check for possible errors
    if(echeck) {
      for(p=0;p<N-1;p++) {
        for(j=0;j<N-1;j++) {
          sums=0;
          for(i=0;i<N-1;i++) {
            sums+=ev(i, p)*tmpm(i, j);
          }
          if(fabs(sums-A_diag(p)*ev(j, p))>ERROR_MARGIN) {
            cout << sums << " " << A_diag(p) << endl;
            cout << "Error, eigenvector check!" << endl;
            return -1;
          }
        }
      }
      if(debug) cout << "ALL OK" << endl;
    }

    double factor=1/((double)N*friction_av);

    if(debug) cout << "Factor is " << factor << endl;
    if(debug) cout << "Calculating c(i,j)" << endl;

    /* This is the old algorithm! I have made it faster...
		for(i=0;i<N;i++) {
			for(j=0;j<N-1;j++) {
				c(i,j)=0;

				for(x=1;x<=i;x++) {
					c(i,j)+=ev(j,x-1);
				}

				for(x=1;x<N;x++) { // k
					for(y=x;y<N;y++) { // l
						c(i,j)-=factor*friction(y)*ev(j,x-1);
					}
				}
			}
		}
     */

    // fast implementation of calculating c(i,j)
    // Die Indizes der Eigenvektoren sind, im Vergleich zum Dokument, vertauscht!
    ap::real_1d_array cij_factor;
    cij_factor.setbounds(0, N-2);
    for(i=0;i<N-1;i++) {
      cij_factor(i)=0.0;
      //cout << "F " << i << endl;
      for(x=1;x<N;x++) { // k
        for(y=x;y<N;y++) { // l
          cij_factor(i)+=factor*friction(y)*ev(x-1,i);
          //cout << ev(i,x-1) << " " << factor*friction(y)*ev(i,x-1) << endl;
        }
      }
    }

    for(i=0;i<N;i++) {
      for(j=0;j<N-1;j++) {
        c(i,j)=0;

        for(x=1;x<=i;x++) {
          c(i,j)+=ev(x-1,j);
        }

        c(i,j)-=cij_factor(j);
      }
    }

    /*
     * DIES IST DIE VERSION VON WISCHI, SIE IST ANDERS ALS DAS DOKUMENT!
		for(i=0;i<N-1;i++) {
			for(j=0;j<N-1;j++) {
				c(i,j)=0;

				for(x=i+1;x<N-1;x++) {
					c(i,j)-=ev(x,j);
				}

				for(x=0;x<N-1;x++) { // k
					for(y=0;y<=x;y++) { // l
						c(i,j)+=factor*friction(y)*ev(x,j);
					}
				}
			}
		}

		cout << "c(i,j)" << endl;
		cout << endl;
		print_matrix(c, N-1);
		cout << endl;
     */

	// debug output
    if(debug) {
      cout << "Printing c(i,j)" << endl;
      cout << endl;
      cout << "N: " << N << endl;
      for(i=0;i<N;i++) {
        for(j=0;j<N-1;j++) {
          cout << " c" << i << "," << j;
        }
        cout << endl;
      }
      cout << "------" << endl;
      for(i=0;i<N-1;i++) {
        cout << " " << cij_factor(i);
      }
      cout << endl;
      cout << "------" << endl;
      for(i=0;i<N;i++) {
        for(j=0;j<N-1;j++) {
          cout << " " << c(i, j);
        }
        cout << endl;
      }
      cout << endl;
    }

    // check sums (Formula 18)
    if(echeck) {
      if(debug) cout << "Identity check" << endl;

      for(x=0;x<N-1;x++) {
        for(y=0;y<N-1;y++) {
          sums=0;
          for(i=0;i<N-1;i++) {
            sums+=ev(i,x)*ev(i,y);
          }
          if(debug) cout << " " << (int)round(sums);
          if(x==y) {
            if(fabs(sums-1.0)>ERROR_MARGIN) {
              cout << "Error identity check, expected 1" << endl;
              return -1;
            }
          }
          else {
            if(fabs(sums)>ERROR_MARGIN) {
              cout << "Error identity check, expected 0" << endl;
              return -1;
            }
          }
        }
        if(debug) cout << endl;
      }
      if(debug) cout << "ALL OK" << endl;
    }

    // DEBUG TODO
    if(false) {
      if(debug) cout << "Formula 48 check" << endl;
      for(p=0;p<N-1;p++) {
        for(x=0;x<N;x++) {
          for(y=0;y<N;y++) {
            if(x>=y) {
              sums=0;
              for(i=y;i<x;i++) {
                sums+=ev(p,i);
              }
              //							cout << c(x,p)-c(y,p) << " " << sums << endl;
              if(fabs(c(x,p)-c(y,p)-sums)>=ERROR_MARGIN) {
                cout << "Error formula 48 check" << endl;
                return -1;
              }
            }
            else {
              sums=0;
              for(i=x;i<y;i++) {
                sums+=ev(p,i);
              }
              //							cout << c(y,p)-c(x,p) << " " << sums << endl;
              if(fabs(c(y,p)-c(x,p)-sums)>=ERROR_MARGIN) {
                cout << "Error formula 48 check" << endl;
                return -1;
              }
            }
          }
        }
      }
      if(debug) cout << "ALL OK" << endl;
    }

    if(debug) cout << "Calculating MSD, Sinc and Scoh" << endl;
    for(i=0;i<MSD_entries;i++) // i as time
    {
      //cout << "time: " << times(i) << endl << endl;
      for(x=0;x<N;x++) // n
      {
        for(y=0;y<N;y++) // m
        {
          phi(x,y)=0.0;

          for(p=0;p<N-1;p++)
          {
            phi(x,y)+=c(x,p)*c(x,p);
            phi(x,y)+=c(y,p)*c(y,p);
            phi(x,y)-=2.0*c(x,p)*c(y,p)*exp(-A_diag(p)*times(i));
          }

          phi(x,y)*=l*l;

          for(j=0;j<numQ;j++) {
            Scoh_woD(j,i)+=exp(-Qval[j] * Qval[j] / 6 * phi(x,y));
          }

          if(x==y)
          {
            msd_woD(i)+=phi(x,y)/(double)(N-1); // /N to average over beads
          }

          if(!povray) phi(x,y)+=6*D*times(i);

          if(x==y)
          {
            msd(i)+=phi(x,y)/(double)(N-1);
            for(j=0;j<numQ;j++) Sinc(j,i)+=exp(-Qval[j] * Qval[j] / 6 * phi(x,y));
          }

          for(j=0;j<numQ;j++) Scoh(j,i)+=exp(-Qval[j] * Qval[j] / 6 * phi(x,y));
        }
        //cout << endl;
      }

      if(debug) {
        bool draw=false;
        if(i==0) { cout << "t=0" << endl; draw=true; }
        if(i==MSD_entries/2) { cout << "t=1/2" << endl; draw=true; }
        if(i==MSD_entries-1) { cout << "t=infty" << endl; draw=true; }
        if(draw) {
          for(x=0;x<N;x++) // n
          {
            for(y=0;y<N;y++) // m
            {
              cout << " " << phi(x,y);
            }
            cout << endl;
          }
        }
      }
      /*
      for(x=0;x<N;x++) // n
      {
        for(y=0;y<N;y++) // m
        {
          cout << " " << phi(x,y);
        }
        cout << endl;
      }
       */

      // POVRAY
      if(povray) {
        int fr=i;

        char buffer[7];
        int errorcode;
        errorcode=snprintf(buffer, 7, "%6.6d", fr);
        if(errorcode!=6) {
          cout << "Errorcode: " << errorcode << endl;
          exit(-1);
        }
        string time_with_zeros=buffer;

        string name="/local/brodeck/povray/ranrouse01/rr"+time_with_zeros+".inc";
        write_povray(phi, N, (double)N*l*l*2, name, "mymat");
        cout << name << endl;
      }
    }

    // Calculate Rouse times
    if(debug) cout << "Calculating Rouse times..." << endl;
    //sort_vector(A_diag, N-1); // should be sorted!!!

    for(i=0;i<N-1;i++) {
      er_av(i)+=A_diag(i);
    }

    if(debug) {
      cout << "DONE with chain!" << endl;
    }
  }

  cout << "pid " << pid << " is ready and waiting for the rest..." << endl;
  MPI_Barrier(MPI_COMM_WORLD);

  if(root) cout << "Gathering!" << endl;

  // get information from all processors and gather
  MPI_Allreduce(MPI_IN_PLACE, er_av.m_Vec, er_av.m_iVecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, msd.m_Vec, msd.m_iVecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, msd_woD.m_Vec, msd_woD.m_iVecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, Sinc.m_Vec, Sinc.m_iVecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, Scoh.m_Vec, Scoh.m_iVecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, Scoh_woD.m_Vec, Scoh_woD.m_iVecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // output, root process only
  if(root) {
    cout << "Saving!" << endl;

    string filename=prefix+"rt";
    ofstream file(filename.c_str(), ios::out);
    file << "#p t_p" << endl;
    for(i=0;i<N-1;i++) {
      file << (i+1) << " " << (double)chains/er_av(i) << endl; // 1/eigenvalue
    }
    file.close();

    filename=prefix+"msd";
    file.open(filename.c_str(), ios::out);
    file << "#t msd" << endl;
    for(i=0;i<MSD_entries;i++) {
      file << times(i) << " " << msd(i)/(double)chains << endl; // 1/eigenvalue
    }
    file.close();

    filename=prefix+"msd_woD";
    file.open(filename.c_str(), ios::out);
    file << "#t msd" << endl;
    for(i=0;i<MSD_entries;i++) {
      file << times(i) << " " << msd_woD(i)/(double)chains << endl; // 1/eigenvalue
    }
    file.close();

    filename=prefix+"Sinc";
    file.open(filename.c_str(), ios::out);
    file << "#t";
    for(i=0;i<numQ;i++) {
      file << " " << Qval[i];
    }
    file << endl;
    for(i=0;i<MSD_entries;i++) {
      file << times(i);
      for(j=0;j<numQ;j++) {
        file << " " << Sinc(j, i)/Sinc(j,0);
      }
      file << endl;
    }
    file.close();

    filename=prefix+"Scoh";
    file.open(filename.c_str(), ios::out);
    file << "#t";
    for(i=0;i<numQ;i++) {
      file << " " << Qval[i];
    }
    file << endl;
    for(i=0;i<MSD_entries;i++) {
      file << times(i);
      for(j=0;j<numQ;j++) {
        file << " " << Scoh(j, i)/Scoh(j,0);
      }
      file << endl;
    }
    file.close();

    filename=prefix+"Scoh_woD";
    file.open(filename.c_str(), ios::out);
    file << "#t";
    for(i=0;i<numQ;i++) {
      file << " " << Qval[i];
    }
    file << endl;
    for(i=0;i<MSD_entries;i++) {
      file << times(i);
      for(j=0;j<numQ;j++) {
        file << " " << Scoh_woD(j, i)/Scoh_woD(j,0);
      }
      file << endl;
    }
    file.close();
  }

  delete[] inputfrict;
  MPI_Barrier(MPI_COMM_WORLD);

  cout << " --> Zeiterfassung Logout: I am processor " << pid << " and I am DONE" << endl << flush;

  if(root) {
    cout << "Current time: " << std::flush;
    system("date");
  }

  if(root) { cout << "MPI_Finalize!" << endl; }

  MPI_Finalize();

  return 0;
}

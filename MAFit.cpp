#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<limits>
#include<regex>
#include<ctime>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define SPLIT_MAX 3      /* in function split */
#define DATANUM_MAX 10000
#define FILENUM_MAX 100

/* in function getM1_ma */
/* this value is max error in M1 calculation */
#define DELTA 1e-8

using namespace std;

/* for fitting */
double rate_n = -1e-2;
double rate_Tr = -1e-2;
double rate_dCp = -1e-2;
double rate_b = -1e-3;

/* effective digit */
int keta_T = 3;
int keta_dH = 6;
int keta_dG = 6;
int keta_M1 = 4;
int keta_dM1dT = 6;
int keta_CpLipid = 6;
int keta_1stTerm = 8;
int keta_2ndTerm = 8;
int keta_3rdTerm = 8;
int keta_CpObs = 6;
int keta_n = 4;
int keta_Tr = 4;
int keta_dCp = 4;
int keta_b = 4;
int keta_J = 10;
int keta_va = 4;
int keta_vb = 4;

/* calculation method : ma(mass action) pps(pseudo phase separation) */
string model = "ma";

string files[FILENUM_MAX];
double M_arr[FILENUM_MAX];
int fileCount = 0;

/* an array of pointers */
/* each pointer points experiment data array */
double** exp_data;

int lengthT;
double* T;
double* celsius;

/* a pointer which points current argument */
int argCount = 1;

/* gas constant (J/K mol) */
double R = 8.31;

/* volume of DSC cell (ml) */
double V_cell = 0.3;

bool move_n = false;
bool move_Tr = false;
bool move_dCp = false;
bool move_b = false;

bool M_specified = false;
bool n_specified = false;
bool Tr_specified = false;
bool dCp_specified = false;
bool b_specified = false;
string Cp1_specified = "none";
string Cp2_specified = "none";

/* parameter range */
double start_n = 10;
double end_n = 30;
double d_n = 0.2;

double start_Tr = 300;
double end_Tr = 330;
double d_Tr = 1;

double start_dCp = 340;
double end_dCp = 360;
double d_dCp = 1;

double start_b = 25;
double end_b = 30;
double d_b = 0.1;

/* extract M value from specified filename */
bool M_auto = false;

/* plot specified file and calculated results */
bool plot = false;

int tmpFileCount = 1;

/* temperature range */
double startCelsius = 5;
double endCelsius = 99;
double dT = 0.05;
double startT = startCelsius + 273.15;
double endT = endCelsius + 273.15;

double M;
double Cp1; /* specific heat of monomer */
double Cp2; /* specific heat of micelle */
string Cp1_file;
string Cp2_file; 
double* Cp1_file_arr;
double* Cp2_file_arr;
double nc = 6;

/* Volume is in unit angstrom^3 */
/* va is temperature dependent fraction and 
 * vb is temperature independent fraction
 * such that v(one molcule) = va*T + vb 
 */
double va = 2 * nc * (1/32.0);   
double vb = 344 + 26 * 2 * nc;

double* Cw_arr;
double* Dw_arr;
double* Vl_arr;
double* Gl_arr;
double* V_arr;
double** Gl_arrs;
double** V_arrs;

void usage();
void waitEnterKey();
string getNextArg(int, char**);
double NaN();
bool isNaN(double);
bool isNum(string);
bool exists(string);
void split(string, string*);
string getTmpFile();
int readFile(string, double*, double*);
double extractConcentration(string);
void integrate(double*, double*);
double getY(double, double*, double*, int);
void getdCp(double, double*);
void getdCp_dCpnum(double, double*);
void getdCp_Cp1num_Cp2num(double*);
void getdCp_Cp1num_Cp2file(double*);
void getdCp_Cp1file_Cp2num(double*);
void getdCp_Cp1file_Cp2file(double*);
void getCp2(double, double*);
void getCp2_Cp2num(double*);
void getCp2_Cp2File(double*);
void getCp2_dCpnum_Cp1num(double, double*);
void getCp2_dCpnum_Cp1file(double, double*);
void readExpData(string, double*);
void getdH(double*, double, double*);
void getdG(double*, double, double*);
void getTds(double*, double*, double*);
void getM1_pps(double*, double, double*);
void getM1_ma(double*, double, double, double*);
void getdM1dT(double*, double*);
void getCpLipid(double, double*, double*, double*, double*, double*, double*);
void get3Terms(double, double*, double*, double*, double*, double*, double*, double*, double*);
double getJ(double*, double*);
double getdJdn(double, double, double*, double, double*);
double getdJdTr(double, double, double*, double, double*);
double getdJddCp(double, double, double*, double, double*);
double getdJdb(double, double, double*, double, double*);
void getCpObs(double*, double*, double*);
double nc2MW(double);
void getCw();
void getDw();
void getVl();
void getGl(double, double*);
void getV(double*, double*);
void printPlotFormat(double, double, double, double);
void printModelFormat(double, double, double, double);
void printJFormat(double, double, double, double);
void printFitFormat(double, double, double, double);
void addConstToArray(double, double*, int);
void copyArray(double*, double*, int);
void printArr(double*);

int main(int argc, char* argv[]) {

  cout.setf(ios::scientific);

  string dsc_cmd = "";
  double n;   /* aggregation number */
  double Tr;  /* temperature at which dH=0 */
  double dCp; /* difference between Cp1 and Cp2 (dCp=Cp1-Cp2) */
  double b;   /* const of integration */

  /* get commandline options */
  
  string arg;
  argCount = 1;
  while (argCount<argc) {
    arg = getNextArg(argc, argv);
    if (arg.substr(0, 1) == "-") {
      if (arg == "-model") {
	arg = getNextArg(argc, argv);
	if (arg == "ma") {
	  model = "ma";
	} else if (arg == "pps") {
	  model = "pps";
	} else {
	  cerr << "Error: Unknown argument " << arg << endl;
	  exit(0);
	}
      } else if (arg == "-M") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  M = stod(arg);
	  if (M <= 0) {
	    cerr << "Error: M must be positive number" << endl;
	    exit(0);
	  }
	  M_specified = true;
	} else if (arg == "auto") {
	  M_auto = true;
        } else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-F") {
	arg = getNextArg(argc, argv);
	files[fileCount] = arg;
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  M_arr[fileCount] = stod(arg);
	  if (M_arr[fileCount] <= 0) {
	    cerr << "Error: M=" << M_arr[fileCount] << " must be positive number " << endl;
	    exit(0);
	  }
	} else {
	  cerr << "Error: Argument " << arg << " does not look line number" << endl;
	  exit(0);
	}
	fileCount++;
      } else if (arg == "-n") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  n = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
	if (n <= 0) {
	  cerr << "Error: n must be positive number" << endl;
	  exit(0);
	}
	n_specified = true;
      } else if (arg == "-Tr") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  Tr = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
	Tr_specified = true;
      } else if (arg == "-dCp") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  dCp = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
	dCp_specified = true;
      } else if (arg == "-b") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  b = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
	b_specified = true;
      } else if (arg == "-Cp1") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  Cp1 = stod(arg);
	  Cp1_specified = "num";
	} else if (exists(arg)) {
	  Cp1_file = arg;
	  Cp1_specified = "file";
	} else {
	  cerr << "Error: File " << arg << " does not exist" << endl;
	  exit(0);
	}
      } else if (arg == "-Cp2") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  Cp2 = stod(arg);
	  Cp2_specified = "num";
	} else if (exists(arg)) {
	  Cp2_file = arg;
	  Cp2_specified = "file";
	} else {
	  cerr << "Error: File " << arg << " does not exist" << endl;
	  exit(0);
	}
      } else if (arg == "-sT") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  startCelsius = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-eT") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  endCelsius = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number" << endl;
	  exit(0);
	}
      } else if (arg == "-dT") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  dT = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look line nuber " << endl;
	  exit(0);
	}
      } else if (arg == "-move") {
	arg = getNextArg(argc, argv);
	if (arg == "n") {
	  move_n = true;
	} else if (arg == "Tr") {
	  move_Tr = true;
	} else if (arg == "dCp") {
	  move_dCp = true;
	} else if (arg == "b") {
	  move_b = true;
	} else {
	  cerr << "Error: Unknown argument " << arg << endl;
	  exit(0);
	}
      } else if (arg == "-sn") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  start_n = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-sTr") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  start_Tr = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-sdCp") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  start_dCp = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-sb") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  start_b = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-en") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  end_n = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-eTr") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  end_Tr = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-edCp") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  end_dCp = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-eb") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  end_b = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-dn") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  d_n = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-dTr") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  d_Tr = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-ddCp") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  d_dCp = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-db") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  d_b = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-nc") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  nc = stod(arg);
	  if (nc <= 0) {
	    cerr << "Error: nc must be positive integer" << endl;
	    exit(0);
	  }
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-va") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  va = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-vb") {
	arg = getNextArg(argc, argv);
	if (isNum(arg)) {
	  vb = stod(arg);
	} else {
	  cerr << "Error: Argument " << arg << " does not look like number " << endl;
	  exit(0);
	}
      } else if (arg == "-plot") {
	plot = true;
      } else if (arg == "-cmd") {
	arg = getNextArg(argc, argv);
	while (arg != "end") {
	  dsc_cmd = dsc_cmd + " " + arg;
	  arg = getNextArg(argc, argv);
	}
      } else if (arg == "-h") {
	usage();
      } else {
	cerr << "Error: Unknown option " << arg << endl;
	exit(0);
      }
    } else {
      files[fileCount] = arg;
      fileCount++;
    }
  }





  
  
  /* =======================preprocess====================== */

  if (startCelsius >= endCelsius) {
    cerr << "Error: start T=" << startCelsius << " is greater than end T=" << endT << endl;
    exit(0);
  }
  if (dT <= 0) {
    cerr << "Error: dT=" << dT << " must be positive number" << endl;
    exit(0);
  }
  if (start_n >= end_n) {
    cerr << "Error: start n=" << start_n << " is greater than end n=" << end_n << endl;
    exit(0);
  }
  if (d_n <= 0) {
    cerr << "Error: dn=" << d_n << " must be positive number" << endl;
    exit(0);
  }
  if (start_Tr >= end_Tr) {
    cerr << "Error: start Tr=" << start_Tr << " is greater than end Tr=" << end_Tr << endl;
    exit(0);
  }
  if (d_Tr <= 0) {
    cerr << "Error: dTr=" << d_Tr << " must be positive number" << endl;
    exit(0);
  }
  if (start_dCp >= end_dCp) {
    cerr << "Error: start dCp=" << start_dCp << " is greater than end dCp=" << end_dCp << endl;
    exit(0);
  }
  if (d_dCp <= 0) {
    cerr << "Error: ddCp=" << d_dCp << " must be positive number" << endl;
    exit(0);
  }
  if (start_b >= end_b) {
    cerr << "Error: start b=" << start_b << " is greater than end b=" << end_b << endl;
    exit(0);
  }
  if (d_b <= 0) {
    cerr << "Error: db=" << d_b << " must be positive number" << endl;
    exit(0);
  }
  if (start_n <= 0) {
    cerr << "Error: start_n=" << start_n << " must be positive number" << endl;
    exit(0);
  }
  if (start_Tr < 0) {
    cerr << "Error: start_Tr=" << start_Tr << " must be positive number" << endl;
    exit(0);
  }

  for (int i=0; i<fileCount; i++) {
    if (!exists(files[i])) {
      cerr << "Error: File " << files[i] << " does not exist" << endl;
      exit(0);
    }
  }


  /* Generates T-array and celsius-array from startT, endT and dT */
  /* Functions in this program (such as dH dG M1 CpLipid etc.) 
   * have same length as T-array and their ith value corresponds to
   * T-array's ith value.
   */
  
  startT = startCelsius + 273.15;
  endT    =  endCelsius + 273.15;
  lengthT = (endCelsius - startCelsius) / dT + 1;
  celsius = new double[lengthT];
  T       = new double[lengthT];
  for (int i=0; i<lengthT; i++) {
    celsius[i] = startCelsius + dT*i;
    T[i]       = startT       + dT*i;
  }


  /* read experiment data files */
  exp_data = new double*[fileCount];
  for (int i=0; i<fileCount; i++) {
    exp_data[i] = new double[lengthT];
    readExpData(files[i], exp_data[i]);
  }


  
  double Cp1_arr[lengthT];
  double Cp2_arr[lengthT];
  double dCp_arr[lengthT];
  double dH_arr[lengthT];
  double dG_arr[lengthT];
  double M1_arr[lengthT];
  double dM1dT_arr[lengthT];
  double CpLipid_arr[lengthT];
  double firstTerm_arr[lengthT];
  double secondTerm_arr[lengthT];
  double thirdTerm_arr[lengthT];
  double CpObs_arr[lengthT];

  Cp1_file_arr = new double[lengthT];
  Cp2_file_arr = new double[lengthT];
  Cw_arr = new double[lengthT];
  Dw_arr = new double[lengthT];
  Vl_arr = new double[lengthT];
  Gl_arr = new double[lengthT];
  V_arr  = new double[lengthT];
  Gl_arrs = new double*[fileCount];
  V_arrs = new double*[fileCount];

  for(int i=0; i<fileCount; i++) {
    Gl_arrs[i] = new double[lengthT];
    V_arrs[i]  = new double[lengthT];
  }


  /* read Cp1 file and Cp2 file */
  double y;
  if (Cp1_specified == "file" && Cp2_specified == "file") {
    double Cp1_X[DATANUM_MAX];
    double Cp1_Y[DATANUM_MAX];
    double Cp2_X[DATANUM_MAX];
    double Cp2_Y[DATANUM_MAX];
    int dataNumCp1 = readFile(Cp1_file, Cp1_X, Cp1_Y);
    int dataNumCp2 = readFile(Cp2_file, Cp2_X, Cp2_Y);
    for (int i=0; i<lengthT; i++) {
      y = getY(celsius[i], Cp1_X, Cp1_Y, dataNumCp1);
      Cp1_file_arr[i] = y;
      y = getY(celsius[i], Cp2_X, Cp2_Y, dataNumCp2);
      Cp2_file_arr[i] = y;
    }
  } else if (Cp1_specified == "file") {
    double Cp1_X[DATANUM_MAX];
    double Cp1_Y[DATANUM_MAX];
    int dataNumCp1 = readFile(Cp1_file, Cp1_X, Cp1_Y);
    for (int i=0; i<lengthT; i++) {
      y = getY(celsius[i], Cp1_X, Cp1_Y, dataNumCp1);
      Cp1_file_arr[i] = y;
    }
  } else if (Cp2_specified == "file") {
    double Cp2_X[DATANUM_MAX];
    double Cp2_Y[DATANUM_MAX];
    int dataNumCp2 = readFile(Cp2_file, Cp2_X, Cp2_Y);
    for (int i=0; i<lengthT; i++) {
      y = getY(celsius[i], Cp2_X, Cp2_Y, dataNumCp2);
      Cp2_file_arr[i] = y;
    }
  }

  /* set concentration */
  if (M_specified) {
    for (int i=0; i<fileCount; i++) {
      M_arr[i] = M;
    }
  } else if (M_auto) {
    for (int i=0; i<fileCount; i++) {
      M_arr[i] = extractConcentration(files[i]);
    }
  }
  
  /* get Cw and Dw */
  getCw();
  getDw();
  getVl();
  getGl(M, Gl_arr);
  getV(V_arr, Gl_arr);
  for (int i=0; i<fileCount; i++) {
    getGl(M_arr[i], Gl_arrs[i]);
    getV(V_arrs[i], Gl_arrs[i]);
  }

  /* print Dw Cw Vl Gl V */
#ifdef TESTV
  cerr << "1:celsius 2:Cw 3:Dw 4:Vl 5:Gl 6:V 7:Cw*Dw*V" << endl;
  cerr << "M=" << M << endl;
  double CwDwV;
  for (int i=0; i<lengthT; i++) {
    CwDwV = Cw_arr[i] * Dw_arr[i] * V_arr[i];
    cout << celsius[i] << "\t";
    cout << Cw_arr[i] << "\t";
    cout << Dw_arr[i] << "\t";
    cout << Vl_arr[i] << "\t";
    cout << Gl_arr[i] << "\t";
    cout << V_arr[i] << "\t";
    cout << CwDwV << "\n";
  }
  exit(0);
#endif
  
  /* measure calulation time */
#ifdef TIME

  n = 20;
  Tr = 315;
  dCp = 340;
  b = 29;
  M = 15;
  Cp1 = 1500;
  getdCp_dCpnum(dCp, dCp_arr);
  getCp2_dCpnum_Cp1num(dCp, Cp2_arr);
  exp_data[0] = new double[lengthT];
  
  double dH_Time = 0;
  double dG_Time = 0;
  double M1_Time = 0;
  double dM1dT_Time = 0;
  double CpLipid_Time = 0;
  double CpObs_Time = 0;
  double J_Time = 0;
  int iteration = 1000;
  clock_t begin;
  clock_t end;

  cout << "=====================================" << endl;
  cout << "iteration=" << iteration << endl;
  begin = clock();
  for (int i=0; i<iteration; i++) {
    getdH(dCp_arr, Tr, dH_arr);
  }
  end = clock();
  dH_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "dH calculation:      elapsed time=" << dH_Time << "(sec)" << endl;

  begin = clock();
  for (int i=0; i<iteration; i++) {
    getdG(dH_arr, b, dG_arr);
  }
  end = clock();
  dG_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "dG calculation:      elapsed time=" << dG_Time << "(sec)" << endl;

  begin = clock();
  for (int i=0; i<iteration; i++) {
    getM1_ma(dG_arr, M, n, M1_arr);
  }
  end = clock();
  M1_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "M1 calculation:      elapsed time=" << M1_Time << "(sec)" << endl;

  begin = clock();
  for (int i=0; i<iteration; i++) {
     getdM1dT(M1_arr, dM1dT_arr);
  }
  end = clock();
  dM1dT_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "dM1/dT calculation:  elapsed time=" << dM1dT_Time << "(sec)" << endl;

  begin = clock();
  for (int i=0; i<iteration; i++) {
    getCpLipid(M, dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
  }
  end = clock();
  CpLipid_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "CpLipid calculation: elapsed time=" << CpLipid_Time << "(sec)" << endl;

  begin = clock();
  for (int i=0; i<iteration; i++) {
    getCpObs(CpLipid_arr, V_arr, CpObs_arr);
  }
  end = clock();
  CpObs_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "CpObs calculation:   elapsed time=" << CpObs_Time << "(sec)" << endl;

  begin = clock();
  for (int i=0; i<iteration; i++) {
      getJ(CpObs_arr, exp_data[0]);
  }
  end = clock();
  J_Time = double(end - begin) / CLOCKS_PER_SEC;
  cout << "J calculation:       elapsed time=" << J_Time << "(sec)" << endl;

  cout << "=====================================" << endl;
  double all_Time = dH_Time + dG_Time + M1_Time + dM1dT_Time + CpLipid_Time + CpObs_Time + J_Time;
  cout << "dH calculation time      " << dH_Time / all_Time * 100 << "(%)" << endl;
  cout << "dG calculation time      " << dG_Time / all_Time * 100 << "(%)" << endl;
  cout << "M1 calculation time      " << M1_Time / all_Time * 100 << "(%)" << endl;
  cout << "dM1dT calculation time   " << dM1dT_Time / all_Time * 100 << "(%)" << endl;
  cout << "CpLipid calculation time " << CpLipid_Time / all_Time * 100 << "(%)" << endl;
  cout << "CpObs calculation time   " << CpObs_Time / all_Time * 100 << "(%)" << endl;
  cout << "J calculation time       " << J_Time / all_Time * 100 << "(%)" << endl;
  
  exit(0);
#endif




  /* ==========================main======================== */
  string object;
  bool move_parameter;
  
  if (!move_n && !move_Tr && !move_dCp && !move_b) {
    move_parameter = false;
  } else {
    move_parameter = true;
  }

  if (fileCount == 0 && move_parameter) {
    cerr << "Error: Specify at least one file in order to move parameter and calculate J" << endl;
    exit(0);
  } else if (plot && fileCount == 0) {
    cerr << "Error: Specify at least one file in order to plot calculated results" << endl;
    exit(0);
  } else if (plot) {
    object = "plot";
  } else if (fileCount == 0 && !move_parameter) {
    object = "calc_model";
  } else if (fileCount > 0 && move_parameter) {
    object = "calc_J";
  } else if (fileCount > 0 && !move_parameter) {
    object = "fitting";
  }

  if (object == "plot") {
    printPlotFormat(n, Tr, dCp, b);
  } else if (object == "calc_model") {
    printModelFormat(n, Tr, dCp, b);
  } else if (object == "calc_J") {
    printJFormat(n, Tr, dCp, b);
  } else if (object == "fitting") {
    printFitFormat(n, Tr, dCp, b);
  }



  if (object == "plot") {
    string tmpFile;
    string tmpFiles = "";
    string plotFiles = "";
    
    for (int j=0; j<fileCount; j++) {
      getCp2(dCp, Cp2_arr);
      getdCp(dCp, dCp_arr);	
      getdH(dCp_arr, Tr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[j], M1_arr);
      } else {
	getM1_ma(dG_arr, M_arr[j], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[j], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      get3Terms(M_arr[j], dCp_arr, dH_arr, M1_arr, dM1dT_arr, Cp2_arr, firstTerm_arr, secondTerm_arr, thirdTerm_arr);
      getCpObs(CpLipid_arr, V_arrs[j], CpObs_arr);

      tmpFile = getTmpFile();
      ofstream out(tmpFile);
      out.setf(ios::scientific);
      out << "#celsius dH     dG     M1     dM1/dT     CpLipid 1stTerm 2ndTerm 3rdTerm CpObs" << endl;
      for (int i=0; i<lengthT; i++) {
	out << setprecision(keta_T) << celsius[i] << "\t";
	out << setprecision(keta_dH) << dH_arr[i] << "\t";
	out << setprecision(keta_dG) << dG_arr[i] << "\t";
	out << setprecision(keta_M1) << M1_arr[i] << "\t";
	out << setprecision(keta_dM1dT) << dM1dT_arr[i] << "\t";
	out << setprecision(keta_CpLipid) << CpLipid_arr[i] << "\t";
	out << setprecision(keta_1stTerm) << firstTerm_arr[i] << "\t";
	out << setprecision(keta_2ndTerm) << secondTerm_arr[i] << "\t";
	out << setprecision(keta_3rdTerm) << thirdTerm_arr[i] << "\t";
	out << setprecision(keta_CpObs) << CpObs_arr[i] << endl;
      }
      out.close();
      tmpFiles = tmpFiles + tmpFile + " ";
      plotFiles = plotFiles + files[j] + " ";
    }
    
    string cmd = "/home/hatatani/bin/dsc_graph2 " + plotFiles + " " + tmpFiles + " -use1-" + to_string(fileCount) + " 1:2 -use" + to_string(fileCount+1) + "-" + to_string(fileCount*2) + " 1:10 -lc1-" + to_string(fileCount) + " 1 -lc" + to_string(fileCount+1) + "-" + to_string(fileCount*2) + " 3 " + dsc_cmd;
    const char* _cmd = cmd.c_str();
    system(_cmd);

    exit(0);
  }
  

  
  if (object == "calc_model") {
    
    getCp2(dCp, Cp2_arr);
    getdCp(dCp, dCp_arr);	
    getdH(dCp_arr, Tr, dH_arr);
    getdG(dH_arr, b, dG_arr);
    if (model == "pps") {
      getM1_pps(dG_arr, M, M1_arr);
    } else {
      getM1_ma(dG_arr, M, n, M1_arr);
    }
    getdM1dT(M1_arr, dM1dT_arr);
    getCpLipid(M, dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
    get3Terms(M, dCp_arr, dH_arr, M1_arr, dM1dT_arr, Cp2_arr, firstTerm_arr, secondTerm_arr, thirdTerm_arr);
    getCpObs(CpLipid_arr, V_arr, CpObs_arr);

    cout << "#celsius dH     dG     M1     dM1/dT     CpLipid 1stTerm 2ndTerm 3rdTerm CpObs" << endl;
    for (int i=0; i<lengthT; i++) {
      cout << setprecision(keta_T) << celsius[i] << "\t";
      cout << setprecision(keta_dH) << dH_arr[i] << "\t";
      cout << setprecision(keta_dG) << dG_arr[i] << "\t";
      cout << setprecision(keta_M1) << M1_arr[i] << "\t";
      cout << setprecision(keta_dM1dT) << dM1dT_arr[i] << "\t";
      cout << setprecision(keta_CpLipid) << CpLipid_arr[i] << "\t";
      cout << setprecision(keta_1stTerm) << firstTerm_arr[i] << "\t";
      cout << setprecision(keta_2ndTerm) << secondTerm_arr[i] << "\t";
      cout << setprecision(keta_3rdTerm) << thirdTerm_arr[i] << "\t";
      cout << setprecision(keta_CpObs) << CpObs_arr[i] << endl;
    }
  
  exit(0);
  }




  /* move parameters and calculate J */
  if (object == "calc_J") {
    
    double n_min;
    double Tr_min;
    double dCp_min;
    double b_min;
    double J_min = 1e+10;
    double J_single;
    double J_global;

    if (!move_n) {
      start_n = n;
      end_n = n;
    }
    if (!move_Tr) {
      start_Tr = Tr;
      end_Tr = Tr;
    }
    if (!move_dCp) {
      start_dCp = dCp;
      end_dCp = dCp;
    }
    if (!move_b) {
      start_b = b;
      end_b = b;
    }
		
    for (n=start_n; n<=end_n; n+=d_n) {
      for (Tr=start_Tr; Tr<=end_Tr; Tr+=d_Tr) {
	for (dCp=start_dCp; dCp<=end_dCp; dCp+=d_dCp) {
	  for (b=start_b; b<=end_b; b+=d_b) {
	    
	    J_single = 0;
	    J_global = 0;
	    
	    for (int i=0; i<fileCount; i++) {
	      getCp2(dCp, Cp2_arr);
	      getdCp(dCp, dCp_arr);
	      getdH(dCp_arr, Tr, dH_arr);
	      getdG(dH_arr, b, dG_arr);
	      getM1_ma(dG_arr, M_arr[i], n, M1_arr);
	      getdM1dT(M1_arr, dM1dT_arr);
	      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
	      getCpObs(CpLipid_arr, V_arrs[i], CpObs_arr);
	      J_single = getJ(CpObs_arr, exp_data[i]);
	      J_global += J_single;
	    }
	    
	    cout << setprecision(keta_n)  << n << "\t";
	    cout << setprecision(keta_Tr)  << Tr << "\t";
	    if (Cp1_specified != "none" && Cp2_specified != "none" && !move_dCp && !dCp_specified) {
	      cout << "---\t";
	    } else {
	      cout << setprecision(keta_dCp)  << dCp << "\t";
	    }
	    cout << setprecision(keta_b)  << b << "\t";
	    cout << setprecision(keta_J)  << J_global << endl;

	    if (J_min > J_global) {
	      n_min = n;
	      Tr_min = Tr;
	      dCp_min = dCp;
	      b_min = b;
	      J_min = J_global;
	    }
	  }
	}
      }
    }

    cerr << "# J min and parameters" << endl;
    cerr << setprecision(keta_n) << "n=" << n_min << "\t";
    cerr << setprecision(keta_Tr) << "Tr=" << Tr_min << "\t";
    if (Cp1_specified != "none" && Cp2_specified != "none" && !move_dCp && !dCp_specified) {
      cerr << "dCp=---\t";
    } else {
      cerr << setprecision(keta_dCp) << "dCp=" << dCp_min << "\t";
    }
    cerr << setprecision(keta_b) << "b=" << b_min << "\t";
    cerr << setprecision(keta_J) << "J(min)=" << J_min << endl;
    exit(0);
  }

  
  /* fitting */

  if (object == "fitting") {
    cerr << "Error: Fitting method is currently under construction\n";
    exit(0);

    /* if Cp1 file and Cp2 file are specified, fix dCp */
    if (model == "pps") {
      ;
    } else if (model == "ma") {

      double dJdn;
      double dJdTr;
      double dJddCp;
      double dJdb;
      double J;
      double J_single;
      double J_global;
      double J_previous = 1e15;
      
      while (1) {

	J_single = 0;
	J_global = 0;
	
	for (int i=0; i<fileCount; i++) {
	  getCp2(dCp, Cp2_arr);
	  getdCp(dCp, dCp_arr);
	  getdH(dCp_arr, Tr, dH_arr);
	  getdG(dH_arr, b, dG_arr);
	  getM1_ma(dG_arr, M_arr[i], n, M1_arr);
	  getdM1dT(M1_arr, dM1dT_arr);
	  getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
	  getCpObs(CpLipid_arr, V_arrs[i], CpObs_arr);
	  J_single = getJ(CpObs_arr, exp_data[i]);
	  J_global += J_single;
	}
	
	dJdn = getdJdn(n, Tr, dCp_arr, b, Cp2_arr);
	dJdTr = getdJdTr(n, Tr, dCp_arr, b, Cp2_arr);
	dJddCp = getdJddCp(n, Tr, dCp_arr, b, Cp2_arr);
	dJdb = getdJdb(n, Tr, dCp_arr, b, Cp2_arr);
	cout << setprecision(10) << n << "\t" ;
	cout << setprecision(10) << Tr << "\t";
	cout << setprecision(10) << dCp << "\t";
	cout << setprecision(10) << b << "\t";
	cout << setprecision(12) << J_global << endl;

	n += dJdn * rate_n;
	Tr += dJdTr * rate_Tr;
	dCp += dJddCp * rate_dCp;
	b += dJddCp * rate_b;

	if (J_global > J_previous) {
	  cerr << "J increased from previpus step. Calculation stopped." << endl;
	  exit(0);
	} else {
	  J_previous = J_global;
	}
	
      }
    }
  }
}





/*=======================functions=======================*/
void usage() {
  cout << "========================usage=======================" << endl;
  cout << " Fit DSC data using mass action model" << endl;
  cout << " compliler: clang++ " << endl;
  cout << " source: " << __FILE__ << endl;
  cout << "========================options=====================" << endl;
  cout << " -model model : model=ma(mass action) or pps(pseudo phase separation)" << endl;
  cout << " -M value : M value (concentration of data)" << endl;
  cout << "            Above option sets all data's concentrations to specified value" << endl;
  cout << "            If you want to set concentration individually, use -F option" << endl;
  cout << " -M auto  : extract M value from specified filename" << endl;
  cout << " -F file concentration : specify data and concentration of the data" << endl; 
  cout << " -n value   : n value (aggregation number)" << endl;
  cout << " -Tr value  : Tr value" << endl;
  cout << " -dCp value : dCp value" << endl;
  cout << " -b value   : b value" << endl;
  cout << " -Cp1 value or file : Cp1 value or Cp1 file" << endl;
  cout << " -Cp2 value or file : Cp2 value or Cp2 file" << endl;
  cout << " -sT T : start temperature" << endl;
  cout << " -eT T : end temperature" << endl;
  cout << " -dT T : temperature interval" << endl;
  cout << " -move parameter : move parameter to calculate J" << endl;
  cout << "                   (parameter=[n/Tr/dCp/b]) " << endl;
  cout << " -sn start   : set start n" << endl;
  cout << " -sTr start  : set start Tr" << endl;
  cout << " -sdCp start : set start dCp" << endl;
  cout << " -sb start   : set start b" << endl;
  cout << " -en end   : set end n " << endl;
  cout << " -eTr end  : set end Tr " << endl;
  cout << " -edCp end : set end dCp " << endl;
  cout << " -eb end   : set end b " << endl;
  cout << " -dn d     : set interval n" << endl;
  cout << " -dTr d    : set interval Tr" << endl;
  cout << " -ddCp d   : set interval dCp" << endl;
  cout << " -db d     : set interval b" << endl;
  cout << " -plot     : plot file and calculated results" << endl;
  cout << " -cmd [cmd] end : enter dsc_graph2 raw command" << endl;
  cout << "                : (put string \"end\" to the tail of the command)" << endl;
  cout << " -h  : show this help" << endl; 
  cout << "***** below parameters are used to calculate CpObs *****" << endl;
  cout << " -nc value : carbon number " << endl;
  cout << " -va value : v(temperature dependent)" << endl;
  cout << " -vb value : v(temperature independent)" << endl;
  cout << "             v(molecule) = va*T + vb" << endl;
  exit(0);
}

void waitEnterKey() {
  if (cin.get()) {
    return ;
  }
}

string getNextArg(int argc, char* argv[]) {
  string arg;
  if (argCount < argc) {
    arg = string(argv[argCount]);
    argCount++;
  } else {
    cerr << "Error: Next arg is void in function getNextArg" << endl;
    exit(0);
  }
  return arg;
}

double NaN() {
  return numeric_limits<double>::quiet_NaN();
}

bool isNaN(double num1) {
  double num2 = num1;
  if (num1 != num2) {
    return true;
  } else {
    return false;
  }
}

bool isNum(string str) {
  try {
    double num = stod(str);
    return 1;
  } catch (...) {
    return 0;
    }
}

/* does the given file exist? */
bool exists(string filename) {
  const char *file = filename.c_str();
  FILE *fp = fopen(file, "r");
  if (fp == NULL) {
    return false;
  } else {
    fclose(fp);
    return true;
  }
}

/* split string by space or tab */
void split(string str, string ret[]) {
  int spos, tpos;
  int count = 0;
  string first, last;
  
  while (count < SPLIT_MAX) {
    spos = str.find_first_of(" ");
    tpos = str.find_first_of("\t");
    if (spos == string::npos && tpos == string::npos) {
      ret[count] = str;
      break;
    }
    
    if (tpos == string::npos || (spos != string::npos && tpos != string::npos && spos < tpos)) {
      first = str.substr(0, spos);
      last = str.substr(spos+1);
      if (first != "") {
	ret[count] = first;
	count++;
      }
      str = last;
    } else if (spos == string::npos || (spos != string::npos && tpos != string::npos && tpos < spos)) {
      first = str.substr(0, tpos);
      last = str.substr(tpos+1);
      if (first != "") {
	ret[count] = first;
	count++;
      }
      str = last;
    }
  }
  
}

string getTmpFile() {
  string ret = "/tmp/tmpFile" + to_string(tmpFileCount);
  tmpFileCount++;
  return ret;
}

/* read column 1 and column 2 of given file */
/* returns datanum */
int readFile(string filename, double X[], double Y[]) {
  ifstream file(filename);
  string line;
  string str_arr[SPLIT_MAX];
  int dataCount = 0;

  if (!exists(filename)) {
    cerr << "Error: File " << filename << "does not exist (in function readFile)" << endl;
    exit(0);
  }
  
  while (getline(file, line)) {
    if (line.substr(0, 1) == "#") {continue;}
    split(line, str_arr);
    X[dataCount] = stod(str_arr[0]);
    Y[dataCount] = stod(str_arr[1]);
    dataCount++;
  }
  
  return dataCount;
}

double extractConcentration(string filename) {
  smatch results;
  regex_search(filename, results, regex("[0-9]+mM"));
  if (isNum(results.str())) {
    return stod(results.str());
  } else {
    cerr << "Error: Cannot extract concentration. filename=" << filename << endl;
    exit(0);
  }
}

void integrate(double Y[], double ret[]) {
  double sum = 0;
  ret[0] = 0;
  for (int i=0; i<lengthT-1; i++) {
    sum += (Y[i+1] + Y[i]) * dT / 2;
    ret[i+1] = sum;
  }
}

double getY(double x, double X[], double Y[], int length) {
  double y;
  if (x < X[0] || x > X[length-1]) {
    cerr << "x=" << x << " is out of range in function getY" << endl;
    return NaN();
  }

  for (int i=1; i<length; i++) {
    if (X[i] >= x) {
      y = Y[i-1] + (x - X[i-1]) * (Y[i] - Y[i-1]) / (X[i] - X[i-1]);
      break;
    }
  }
  
  return y;
}

void getdCp(double dCp, double dCp_arr[]) {
  if (move_dCp || dCp_specified) {
    getdCp_dCpnum(dCp, dCp_arr);
  } else if (Cp1_specified == "num" && Cp2_specified == "num") {
    getdCp_Cp1num_Cp2num(dCp_arr);
  } else if (Cp1_specified == "num" && Cp2_specified == "file") {
    getdCp_Cp1num_Cp2file(dCp_arr);
  } else if (Cp1_specified == "file" && Cp2_specified == "num") {
    getdCp_Cp1file_Cp2num(dCp_arr);
  } else if (Cp1_specified == "file" && Cp2_specified == "file") {
    getdCp_Cp1file_Cp2file(dCp_arr);
  }
}

void getdCp_dCpnum(double dCp, double dCp_arr[]) {
  for (int i=0; i<lengthT; i++) {
    dCp_arr[i] = dCp;
  }
}

void getdCp_Cp1num_Cp2num (double dCp_arr[]) {
  for (int i=0; i<lengthT; i++) {
    dCp_arr[i] = Cp1 - Cp2;
  }
}

void getdCp_Cp1num_Cp2file (double dCp_arr[]) {
  for (int i=0; i<lengthT; i++) {
    dCp_arr[i] = Cp1 - Cp2_file_arr[i];
  }
}

void getdCp_Cp1file_Cp2num (double dCp_arr[]) {
  for (int i=0; i<lengthT; i++) {
    dCp_arr[i] = Cp1_file_arr[i] - Cp2;
  }
}

void getdCp_Cp1file_Cp2file (double dCp_arr[]) {
  for (int i=0; i<lengthT; i++) {
    dCp_arr[i] = Cp1_file_arr[i] - Cp2_file_arr[i];
  }
}


void getCp2(double dCp, double Cp2_arr[]) {
  if ((move_dCp || dCp_specified) && Cp1_specified == "num") {
    getCp2_dCpnum_Cp1num(dCp, Cp2_arr);
  } else if ((move_dCp || dCp_specified) && Cp1_specified == "file") {
    getCp2_dCpnum_Cp1file(dCp, Cp2_arr);
  } else if (Cp2_specified == "num") {
    getCp2_Cp2num(Cp2_arr);  
  } else if (Cp2_specified == "file") {
    getCp2_Cp2File(Cp2_arr);
  }
}

void getCp2_Cp2num(double Cp2_arr[]) {
  for (int i=0; i<lengthT; i++) {
    Cp2_arr[i] = Cp2;
  }
}

void getCp2_Cp2File(double Cp2_arr[]) {
  for (int i=0; i<lengthT; i++) {
    Cp2_arr[i] = Cp2_file_arr[i];
  }
}

void getCp2_dCpnum_Cp1num (double dCp, double Cp2_arr[]) {
  for (int i=0; i<lengthT; i++) {
    Cp2_arr[i] = Cp1 - dCp;
  }
}

void getCp2_dCpnum_Cp1file (double dCp, double Cp2_arr[]) {
  for (int i=0; i<lengthT; i++) {
    Cp2_arr[i] = Cp1_file_arr[i] - dCp;
  }
}

void readExpData(string filename, double expData[]) {
  double y;
  double expData_X[DATANUM_MAX];
  double expData_Y[DATANUM_MAX];
  int dataNum = readFile(filename, expData_X, expData_Y);
  for (int i=0; i<lengthT; i++) {
    y = getY(celsius[i], expData_X, expData_Y, dataNum);
    if (isNaN(y)) {
      cerr << "Error: x=" << celsius[i] << " is out of range in function readExpData" << endl;
      exit(0);
    } else {
      expData[i] = y;
    }
  }
}

void getdH(double dCp_arr[], double Tr, double dH_arr[]) {
  double integrated[lengthT];
  integrate(dCp_arr, integrated);
  double y = getY(Tr, T, integrated, lengthT);
  if (isNaN(y)) {
    cerr << "Error: Tr=" << Tr << " is out of range in function getdH" << endl;
    exit(0);
  }
  for(int i=0; i<lengthT; i++) {
    dH_arr[i] = integrated[i] - y;
  }
}

void getdG(double dH_arr[], double b, double dG_arr[]) {
  double integrated[lengthT];
  double dH_T2[lengthT];
  
  for (int i=0; i<lengthT; i++) {
    dH_T2[i] = -1 * dH_arr[i] / (T[i]*T[i]);
  }
  
  integrate(dH_T2, integrated);
  
  for (int i=0; i<lengthT; i++) {
    dG_arr[i] = integrated[i]*T[i] + b*T[i];
  }
}

void getTds(double dH_arr[], double dG_arr[], double Tds_arr[]) {
  for (int i=0; i<lengthT; i++) {
    Tds_arr[i] = dH_arr[i] - dG_arr[i];
  }
}

void getM1_pps(double dG_arr[], double M, double M1_arr[]) {
  for (int i=0; i<lengthT; i++) {
    M1_arr[i] = exp(-1*dG_arr[i] / (R*T[i])) * 1000;
    if (M1_arr[i] > M) {
      M1_arr[i] = M;
    }
  }
}

void getM1_ma(double dG_arr[], double M, double n, double M1_arr[]) {
  double min;
  double max;
  double fx;
  double k;
  
  for (int i=0; i<lengthT; i++) {
    min = 0;
    max = M/1000;
    M1_arr[i] = M/2000;
    do {
      k = pow((1/n)*(M/1000-M1_arr[i]), (1/n)) / (M1_arr[i]);
      if (k == 0) {
	fx = -1;
      } else {
	fx = R*T[i]*log(k) - dG_arr[i];
      }

      if (fx < 0) {
	max = M1_arr[i];
      } else {
	min = M1_arr[i];
      }

      M1_arr[i] = (min+max) / 2;
    } while ((max - min) > (DELTA / 1000));
    M1_arr[i] *= 1000;
  }
}

/*
void getM1_ma(double dG[], double M, double n, double M1[]) {
  double min;
  double max;
  double fx;
  double k;
  
  for (int i=0; i<lengthT; i++) {
    min = 0;
    max = M;
    M1[i] = M / 2;
    do {
      k = pow((1/n)*(M-M1[i]), (1/n)) / M1[i];
      if (k == 0) {
	fx = -1;
      } else {
	fx = R*T[i]*log(k) - dG[i];
      }

      if (fx < 0) {
	max = M1[i];
      } else {
	min = M1[i];
      }

      M1[i] = (min+max) / 2;
    } while ((max - min) > DELTA);
  }
}
*/

void getdM1dT(double M1_arr[], double dM1dT_arr[]) {
  dM1dT_arr[0] = (M1_arr[1] - M1_arr[0]) / dT;
  for (int i=0; i<lengthT-1; i++) {
    dM1dT_arr[i+1] = (M1_arr[i+1] - M1_arr[i]) / dT;
  }
}

/* get heat capacity of lipid (water contribution not included) 
 * (CpLipid) + (water contribution) = (CpObs)
 */
void getCpLipid(double M, double dCp_arr[], double dH_arr[], double Cp2_arr[], double M1_arr[], double dM1dT_arr[], double CpLipid_arr[]) {
  for (int i=0; i<lengthT; i++) {
    CpLipid_arr[i] = 1e-6 * V_cell * (M*Cp2_arr[i] + M1_arr[i]*dCp_arr[i] + dH_arr[i]*dM1dT_arr[i]);
  }
}

/* get 3 components of CpLipid
 * 1stTerm = Cp1*M1   (monomer contribution) 
 * 2ndTerm = Cp2*M2   (micelle contribution)
 * 3rdTerm = dH*dM1dT (monomer <-> micelle contribution)
 */
void get3Terms(double M, double dCp_arr[], double dH_arr[], double M1_arr[], double dM1dT_arr[], double Cp2_arr[], double firstTerm[], double secondTerm[], double thirdTerm[]) {
  double Cp1;
  double M2;
  for (int i=0; i<lengthT; i++) {
    Cp1 = Cp2_arr[i] + dCp_arr[i];
    M2 = M - M1_arr[i];
    firstTerm[i] = 1e-6 * V_cell * M1_arr[i] * Cp1;
    secondTerm[i] = 1e-6 * V_cell * M2 * Cp2_arr[i];
    thirdTerm[i] = 1e-6 * V_cell * dH_arr[i] * dM1dT_arr[i];
  }
}

double getJ(double Cp[], double expData[]) {
  double J = 0;
  for (int i=0; i<lengthT; i++) {
    J += (Cp[i] - expData[i]) * (Cp[i] - expData[i]);
  }
  return J;
}

double getdJdn(double n, double Tr, double dCp_arr[], double b, double Cp2_arr[]) {

  if (model == "pps") {
    cerr << "Error: In pps model, n is not a parameter(in function getdJdn)" << endl;
    exit(0);
  }
  
  double dn = 1e-5;
  double J_single;
  double J_global;
  double J_single_tmp1;
  double J_single_tmp2;
  double J_global_tmp1;
  double J_global_tmp2;
  double dH_arr[lengthT];
  double dG_arr[lengthT];
  double M1_arr[lengthT];
  double dM1dT_arr[lengthT];
  double CpLipid_arr[lengthT];

  getdH(dCp_arr, Tr, dH_arr);
  getdG(dH_arr, b, dG_arr);

  do {
    J_single = 0;
    J_global = 0;
    J_single_tmp1 = 0;
    J_single_tmp2 = 0;
    J_global_tmp1 = 0;
    J_global_tmp2 = 0;
    
    for (int i=0; i<fileCount; i++) {	
      getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single = getJ(CpLipid_arr, exp_data[i]);
      J_global += J_single;
    }
        
    for (int i=0; i<fileCount; i++) {	
      getM1_ma(dG_arr, M_arr[i], n+dn, M1_arr);
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp1 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp1 += J_single_tmp1;
    }

    for (int i=0; i<fileCount; i++) {	
      getM1_ma(dG_arr, M_arr[i], n-dn, M1_arr);
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp2 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp2 += J_single_tmp2;
    }

    dn /= 2;
  } while ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) > 0);

  /*
  if ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) == 0) {
    cerr << "Error: Partial derivative dJ/dn is 0 in function getdJdn. Calculation stopped." << endl;
    exit(0);
  }
   */

  return (J_global_tmp1 - J_global) / dn;
}

double getdJdTr(double n, double Tr, double dCp_arr[], double b, double Cp2_arr[]) {
  
  double dTr = 1e-5;
  double J_single;
  double J_global;
  double J_single_tmp1;
  double J_single_tmp2;
  double J_global_tmp1;
  double J_global_tmp2;
  double dH_arr[lengthT];
  double dG_arr[lengthT];
  double M1_arr[lengthT];
  double dM1dT_arr[lengthT];
  double CpLipid_arr[lengthT];

  do {
    J_single = 0;
    J_global = 0;
    J_single_tmp1 = 0;
    J_single_tmp2 = 0;
    J_global_tmp1 = 0;
    J_global_tmp2 = 0;
    
    for (int i=0; i<fileCount; i++) {	
      getdH(dCp_arr, Tr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single = getJ(CpLipid_arr, exp_data[i]);
      J_global += J_single;
    }
        
    for (int i=0; i<fileCount; i++) {	
      getdH(dCp_arr, Tr+dTr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp1 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp1 += J_single_tmp1;
    }

    for (int i=0; i<fileCount; i++) {	
      getdH(dCp_arr, Tr-dTr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp2 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp2 += J_single_tmp2;
    }

    dTr /= 2;
  } while ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) > 0);

  /*
  if ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) == 0) {
    cerr << "Error: Partial derivative dJ/dTr is 0 in function getdJdTr. Calculation stopped." << endl;
    exit(0);
  }
   */

  return (J_global_tmp1 - J_global) / dTr;
}

double getdJddCp(double n, double Tr, double dCp_arr[], double b, double Cp2_arr[]) {

  if (model == "pps") {
    cerr << "Error: In pps model, n is not a parameter(in function getdJdn)" << endl;
    exit(0);
  }
  
  double ddCp = 1e-5;
  double J_single;
  double J_global;
  double J_single_tmp1;
  double J_single_tmp2;
  double J_global_tmp1;
  double J_global_tmp2;
  double dCp_arr_tmp1[lengthT];
  double dCp_arr_tmp2[lengthT];
  double dH_arr[lengthT];
  double dG_arr[lengthT];
  double M1_arr[lengthT];
  double dM1dT_arr[lengthT];
  double CpLipid_arr[lengthT];

  do {
    J_single = 0;
    J_global = 0;
    J_single_tmp1 = 0;
    J_single_tmp2 = 0;
    J_global_tmp1 = 0;
    J_global_tmp2 = 0;
    copyArray(dCp_arr, dCp_arr_tmp1, lengthT);
    copyArray(dCp_arr, dCp_arr_tmp2, lengthT);
    addConstToArray(ddCp, dCp_arr_tmp1, lengthT);
    addConstToArray(-1*ddCp, dCp_arr_tmp2, lengthT);
    
    for (int i=0; i<fileCount; i++) {	
      getdH(dCp_arr, Tr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single = getJ(CpLipid_arr, exp_data[i]);
      J_global += J_single;
    }
        
    for (int i=0; i<fileCount; i++) {	
      getdH(dCp_arr_tmp1, Tr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp1 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp1 += J_single_tmp1;
    }

    for (int i=0; i<fileCount; i++) {	
      getdH(dCp_arr_tmp2, Tr, dH_arr);
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp2 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp2 += J_single_tmp2;
    }

    ddCp /= 2;
  } while ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) > 0);

  /*
  if ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) == 0) {
    cerr << "Error: Partial derivative dJ/ddCp is 0 in function getdJddCp. Calculation stopped." << endl;
    exit(0);
  }
   */

  return (J_global_tmp1 - J_global) / ddCp;
}

double getdJdb(double n, double Tr, double dCp_arr[], double b, double Cp2_arr[]) {
  
  double db = 1e-5;
  double J_single;
  double J_global;
  double J_single_tmp1;
  double J_single_tmp2;
  double J_global_tmp1;
  double J_global_tmp2;
  double dH_arr[lengthT];
  double dG_arr[lengthT];
  double M1_arr[lengthT];
  double dM1dT_arr[lengthT];
  double CpLipid_arr[lengthT];

  getdH(dCp_arr, Tr, dH_arr);
	
  do {
    J_single = 0;
    J_global = 0;
    J_single_tmp1 = 0;
    J_single_tmp2 = 0;
    J_global_tmp1 = 0;
    J_global_tmp2 = 0;
    
    for (int i=0; i<fileCount; i++) {	
      getdG(dH_arr, b, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single = getJ(CpLipid_arr, exp_data[i]);
      J_global += J_single;
    }
        
    for (int i=0; i<fileCount; i++) {	
      getdG(dH_arr, b+db, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp1 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp1 += J_single_tmp1;
    }

    for (int i=0; i<fileCount; i++) {	
      getdG(dH_arr, b-db, dG_arr);
      if (model == "pps") {
	getM1_pps(dG_arr, M_arr[i], M1_arr);
      } else if (model == "ma") {
	getM1_ma(dG_arr, M_arr[i], n, M1_arr);
      }
      getdM1dT(M1_arr, dM1dT_arr);
      getCpLipid(M_arr[i], dCp_arr, dH_arr, Cp2_arr, M1_arr, dM1dT_arr, CpLipid_arr);
      
      J_single_tmp2 = getJ(CpLipid_arr, exp_data[i]);
      J_global_tmp2 += J_single_tmp2;
    }

    db /= 2;
  } while ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) > 0);

  /*
  if ((J_global_tmp1 - J_global) * (J_global_tmp2 - J_global) == 0) {
    cerr << "Error: Partial derivative dJ/db is 0 in function getdJdb. Calculation stopped." << endl;
    exit(0);
  }
   */

  return (J_global_tmp1 - J_global) / db;
}

/* get DSC thermogram (the one after subtracted water/water baseline) */
void getCpObs(double CpLipid_arr[], double V_arr[], double CpObs_arr[]) {
  for (int i=0; i<lengthT; i++) {
    CpObs_arr[i] = CpLipid_arr[i] - Cw_arr[i]*Dw_arr[i]*V_arr[i];
  }
}

/* get molecular weight from number of carbon */
/* assumes PC head, 2 hydrocarbon chains */
double nc2MW() {
  return 285.24 + 28.05*nc;
}

/* specific heat of water (J/K g) */
void getCw() {
  double T;
  for (int i=0; i<lengthT; i++) {
    T = celsius[i];
    Cw_arr[i] = -1.06185e-07*T*T*T + 2.92334e-05*T*T - 0.00178079*T + 4.209;
  }
}

/* density of water (g/ml) */
void getDw() {
  double T;
  for (int i=0; i<lengthT; i++) {
    T = celsius[i];
    Dw_arr[i] = -3.60686e-06*T*T - 7.07451e-05*T + 1.00069;
  }
}

/* get specific volume of lipid (ml/g) */
void getVl() {
  double v_molecule;
  double MW = nc2MW();
  for (int i=0; i<lengthT; i++) {
    v_molecule = va * celsius[i] + vb;
    Vl_arr[i] = v_molecule * 6.02 / 10 / MW;
  }
}

/* get mass of lipid */
void getGl(double M, double Gl_arr[]) {
  double MW = nc2MW();
  for (int i=0; i<lengthT; i++) {
    Gl_arr[i] = V_cell * M * MW / 1e6;
  }
}

/* get volume of lipid */
void getV(double V_arr[], double Gl_arr[]) {
  for (int i=0; i<lengthT; i++) {
    V_arr[i] = Vl_arr[i] * Gl_arr[i];
  }
}

void printPlotFormat(double n, double Tr, double dCp, double b) {

  if (n_specified && model == "pps") {
    cerr << "Error: In pps model, n is not a parameter" << endl;
    exit(0);
  }
  
  cerr << "==========================================" << endl;
  cerr << " Plot specified experimental data along with calculated results" << endl;
  if (model == "pps") {
    cerr << " Model: PPS model" << endl;
  } else if (model == "ma") {
    cerr << " Model: MA model " << endl;
  }
  cerr << " Specified files and concentration: " << endl;
  for (int i=0; i<fileCount; i++) {
    if (M_arr[i] == 0) {
      cerr << "     " << files[i] << ": UNDEFINED" << endl;
    } else {
      cerr << "     " << files[i] << ": " << M_arr[i] << "mM" << endl;
    }
  }
  cerr << " Temperature range " << startCelsius << "-" << endCelsius << " dT=" << dT << endl;
  cerr << " Output: 1:celsius 2:dH 3:dG 4:M1 5:dM1dT 6:CpLipid 7:1stTerm 8:2ndTerm 9:3rdTerm 10:CpObs" << endl;
  cerr << " 1stTerm = Cp1*M1   2ndTerm = Cp2*M2   3rdTerm = dH*dM1/dT" << endl;
  cerr << "===========Specified parameter============" << endl;
  if (model == "pps") {
    cerr << " n: not used in this model" << endl;
  } else if (n_specified) {
    cerr << " n: " << n << endl;
  } else {
    cerr << " n: undefined" << endl;
  }
  
  if (Tr_specified) {
    cerr << " Tr: " << Tr << endl;
  } else {
    cerr << " Tr: undefined" << endl;
  }
  
  if (dCp_specified) {
    cerr << " dCp: " << dCp << endl;
  } else {
    cerr << " dCp: undefined " << endl;
  }
  
  if (b_specified) {
    cerr << " b: " << b << endl;
  } else {
    cerr << " b: undeifned" << endl;
  }

  
  if (Cp1_specified == "file") {
    cerr << " Cp1: " << Cp1_file << endl;
  } else if (Cp1_specified == "num") {
    cerr << " Cp1: " << Cp1 << endl;
  } else {
    cerr << " Cp1: undefined" << endl;
  }
  
  if (Cp2_specified == "file") {
    cerr << " Cp2: " << Cp2_file << endl;
  } else if (Cp2_specified == "num") {
    cerr << " Cp2: " << Cp2 << endl;
  } else {
    cerr << " Cp2: undefined" << endl;
  }

  cerr << " va: " << va << endl;
  cerr << " vb: " << vb << endl; 

  if (!n_specified && model == "ma") {
    cerr << "Error: Specify n" << endl;
    exit(0);
  } else if (!Tr_specified) {
    cerr << "Error: Specify Tr" << endl;
    exit(0);
  } else if (!b_specified) {
    cerr << "Error: Specify b" << endl;
    exit(0);
  }
  
  int count = 0;
  if (dCp_specified) {count++;}
  if (Cp1_specified != "none") {count++;}
  if (Cp2_specified != "none") {count++;}
  if (count < 2) {
    cerr << "Error: Specify at least two of these: dCp, Cp1 and Cp2" << endl;
    exit(0);
  }

  for (int i=0; i<fileCount; i++) {
    if (M_arr[i] == 0) {
      cerr << "Error: specify concentration of data " << files[i] << endl;
      exit(0);
    }
  }
  
  cerr << "==========================================" << endl;
  cerr << " Press Enter key to start" << endl;
  waitEnterKey();
}

void printModelFormat(double n, double Tr, double dCp, double b) {

  if (n_specified && model == "pps") {
    cerr << "Error: In pps model, n is not a parameter" << endl;
    exit(0);
  }
  
  cerr << "==========================================" << endl;
  cerr << " Calculate dH dG M1 dM1/dT CpLipid 1stTerm 2ndTerm 3rdTerm CpObs" << endl;
  if (model == "pps") {
    cerr << " Model: PPS model" << endl;
  } else if (model == "ma") {
    cerr << " Model: MA model " << endl;
  }
  cerr << " Temperature range " << startCelsius << "-" << endCelsius << " dT=" << dT << endl;
  cerr << " Output: 1:celsius 2:dH 3:dG 4:M1 5:dM1dT 6:CpLipid 7:1stTerm 8:2ndTerm 9:3rdTerm 10:CpObs" << endl;
  cerr << " 1stTerm = Cp1*M1   2ndTerm = Cp2*M2   3rdTerm = dH*dM1/dT" << endl;
  cerr << "===========Specified parameter============" << endl;
  if (model == "pps") {
    cerr << " n: not used in this model" << endl;
  } else if (n_specified) {
    cerr << " n: " << n << endl;
  } else {
    cerr << " n: undefined" << endl;
  }
  
  if (Tr_specified) {
    cerr << " Tr: " << Tr << endl;
  } else {
    cerr << " Tr: undefined" << endl;
  }
  
  if (dCp_specified) {
    cerr << " dCp: " << dCp << endl;
  } else {
    cerr << " dCp: undefined " << endl;
  }
  
  if (b_specified) {
    cerr << " b: " << b << endl;
  } else {
    cerr << " b: undeifned" << endl;
  }

  if (M_specified) {
    cerr << " M: " << M << endl;
  } else {
    cerr << " M: undefined" << endl;
  }
  
  if (Cp1_specified == "file") {
    cerr << " Cp1: " << Cp1_file << endl;
  } else if (Cp1_specified == "num") {
    cerr << " Cp1: " << Cp1 << endl;
  } else {
    cerr << " Cp1: undefined" << endl;
  }
  
  if (Cp2_specified == "file") {
    cerr << " Cp2: " << Cp2_file << endl;
  } else if (Cp2_specified == "num") {
    cerr << " Cp2: " << Cp2 << endl;
  } else {
    cerr << " Cp2: undefined" << endl;
  }

  cerr << " va: " << va << endl;
  cerr << " vb: " << vb << endl; 

  if (!n_specified && model == "ma") {
    cerr << "Error: Specify n" << endl;
    exit(0);
  } else if (!Tr_specified) {
    cerr << "Error: Specify Tr" << endl;
    exit(0);
  } else if (!b_specified) {
    cerr << "Error: Specify b" << endl;
    exit(0);
  } else if (!M_specified) {
    cerr << "Error: Specify M" << endl;
    exit(0);
  }
  
  int count = 0;
  if (dCp_specified) {count++;}
  if (Cp1_specified != "none") {count++;}
  if (Cp2_specified != "none") {count++;}
  if (count < 2) {
    cerr << "Error: Specify at least two of these: dCp, Cp1 and Cp2" << endl;
    exit(0);
  }

  cerr << "==========================================" << endl;
  cerr << " Press Enter key to start" << endl;
  waitEnterKey();
}

void printJFormat(double n, double Tr, double dCp, double b)  {

  if (model == "pps" && (move_n || n_specified)) {
    cerr << "Error: In pps model, n is not a parameter" << endl;
    exit(0);
  }
  
      cerr << "==========================================" << endl;
      cerr << " Calculate J" << endl;
      if (model == "pps") {
	cerr << " Model: PPS mdoel" << endl;
      } else if (model == "ma") {
	cerr << " Model: MA model " << endl;
      }
      cerr << " Move parameter: " << endl;
      if (move_n) {
	cerr << "        n  range " << start_n << " to " << end_n << " d=" << d_n << endl;
      }
      if (move_Tr) {
	cerr << "        Tr range " << start_Tr << " to " << end_Tr << " d=" << d_Tr << endl;
      }
      if (move_dCp) {
	cerr << "        dCp range " << start_dCp << " to " << end_dCp << " d=" << d_dCp << endl;
      }
      if (move_b) {
	cerr << "        b range " << start_b << " to " << end_b << " d=" << d_b << endl;
      }
      cerr << " Specified files and concentration: " << endl;
      for (int i=0; i<fileCount; i++) {
	if (M_arr[i] == 0) {
	  cerr << "     " << files[i] << ": UNDEFINED" << endl;
	} else {
	  cerr << "     " << files[i] << ": " << M_arr[i] << "mM" << endl;
	}
      }
      cerr << " Temperature range " << startCelsius << "-" << endCelsius << " dT=" << dT << endl;
      if (model == "pps") {
	cerr << " Output: 1:Tr 2:dCp 3:b 4:J" << endl;
      } else if (model == "ma") {
	cerr << " Output: 1:n 2:Tr 3:dCp 4:b 5:J" << endl;
      }
      cerr << "===========Specified parameter============" << endl;

      if (model == "pps") {
	cerr << " n: not used in this model" << endl;
      } else if (move_n) {
	cerr << " n: move" << endl;
      } else if (n_specified) {
	cerr << " n: " << n << endl;
      } else {
	cerr << " n: udefined " << endl;
      }

      if (move_Tr) {
	cerr << " Tr: move" << endl;
      } else if (Tr_specified) {
	cerr << " Tr: " << Tr << endl;
      } else {
	cerr << " Tr: undefined" << endl;
      }

      if (move_dCp) {
	cerr << " dCp: move" << endl;
      } else if (dCp_specified) {
	cerr << " dCp: " << dCp << endl;
      } else {
	cerr << " dCp: undefined" << endl;
      }

      if (move_b) {
	cerr << " b: move" << endl;
      } else if (b_specified) {
	cerr << " b: " << b << endl;
      } else {
	cerr << " b: undefined" << endl;
      }

      cerr << " va: " << va << endl;
      cerr << " vb: " << vb << endl;

      if (Cp1_specified == "file") {
	cerr << " Cp1: " << Cp1_file << endl;
      } else if (Cp1_specified == "num") {
	cerr << " Cp1: " << Cp1 << endl;
      } else {
	cerr << " Cp1: undefined" << endl;
      }
    
      if (Cp2_specified == "file" && move_dCp) {
	cerr << " Cp2: " << Cp2_file << " (this file is going to be ignored)" << endl;
      } else if (Cp2_specified == "file") {
      cerr << " Cp2: " << Cp2_file << endl;
      } else if (Cp2_specified == "num" && move_dCp) {
	cerr << " Cp2: " << Cp2 << " (this value is going to be ignored)" << endl;
      } else if (Cp2_specified == "num") {
	cerr << " Cp2: " << Cp2 << endl;
      } else {
	cerr << " Cp2: undefined" << endl;
      }

      if (model == "ma" && !move_n && !n_specified) {
	cerr << "Error: Specify n" << endl;
	exit(0);
      } else if (!move_Tr && !Tr_specified) {
	cerr << "Error: Specify Tr" << endl;
	exit(0);
      } else if (!move_b && !b_specified) {
	cerr << "Error: Specify b" << endl;
	exit(0);
      }

      if (move_dCp) {
	if (Cp1_specified == "none") {
	  cerr << "Error: Specify Cp1" << endl;
	  exit(0);
	}
      } else {
	int count = 0;
	if (dCp_specified) {count++;}
	if (Cp1_specified != "none") {count++;}
	if (Cp2_specified != "none") {count++;}
	if (count < 2) {
	  cerr << "Error: Specify at least two of these: dCp, Cp1 and Cp2" << endl;
	  exit(0);
	}
      }

      for (int i=0; i<fileCount; i++) {
	if (M_arr[i] == 0) {
	  cerr << "Error: specify concentration of data " << files[i] << endl;
	  exit(0);
	}
      }
      
      cerr << "==========================================" << endl;
      cerr << " Press Enter key to start" << endl;
      waitEnterKey();
      
}

void printFitFormat(double n, double Tr, double dCp, double b) {

  if (model == "pps" && n_specified) {
    cerr << "Error: In pps model n is not a parameter" << endl;
    exit(0);
  }
  
    cerr << "==========================================" << endl;
    cerr << " Fitting" << endl;
    if (model == "pps") {
      cerr << " Model: PPS mdoel" << endl;
    } else if (model == "ma") {
      cerr << " Model: MA model " << endl;
    }
    cerr << " Specified files and concentration: " << endl;
    for (int i=0; i<fileCount; i++) {
      if (M_arr[i] == 0) {
	cerr << "     " << files[i] << ": UNDEFINED" << endl;
      } else {
	cerr << "     " << files[i] << ": " << M_arr[i] << "mM" << endl;
      }
    }
    cerr << " Temperature range " << startCelsius << "-" << endCelsius << " dT=" << dT << endl;

    if (model == "pps") {
      cerr << " Output: 1:Tr 2:dCp 3:b 4:J" << endl;
    } else if (model == "ma") {
      cerr << " Output: 1:n 2:Tr 3:dCp 4:b 5:J" << endl;
    }

    cerr << "===========Specified parameter============" << endl;
    if (model == "pps") {
      cerr << " n: not used in this model" << endl;
    } else if (n_specified) {
      cerr << " n: " << n << endl;
    } else {
      cerr << " n: undefined" << endl;
    }

    if (Tr_specified) {
      cerr << " Tr: " << Tr << endl;
    } else {
      cerr << " Tr: undefined" << endl;
    }

    if (dCp_specified) {
      cerr << " dCp: " << dCp << endl;
    } else {
      cerr << " dCp: undefined" << endl;
    }

    if (b_specified) {
      cerr << " b: " << b << endl;
    } else {
      cerr << " b: undefined" << endl;
    }

    cerr << " va: " << va << endl;
    cerr << " vb: " << vb << endl;

    if (Cp1_specified == "file") {
      cerr << " Cp1: " << Cp1_file << endl;
    } else if (Cp1_specified == "num") {
      cerr << " Cp1: " << Cp1 << endl;
    } else {
      cerr << " Cp1: undefined" << endl;
    }
    
    if (Cp2_specified == "file") {
      cerr << " Cp2: " << Cp2_file << endl;
    } else if (Cp1_specified == "num") {
      cerr << " Cp2: " << Cp2 << endl;
    } else {
      cerr << " Cp2: undefined" << endl;
    }

    if (model == "ma" && !n_specified) {
      cerr << "Error: Specify n" << endl;
      exit(0);
    } else if (!Tr_specified) {
      cerr << "Error: Specify Tr" << endl;
      exit(0);
    } else if (!b_specified) {
      cerr << "Error: Specify b" << endl;
      exit(0);
    }

 
    int count = 0;
    if (dCp_specified) {count++;}
    if (Cp1_specified != "none") {count++;}
    if (Cp2_specified != "none") {count++;}
    if (count < 2) {
      cerr << "Error: Specify at least two of these: dCp, Cp1 and Cp2" << endl;
      exit(0);
    }

    for (int i=0; i<fileCount; i++) {
      if (M_arr[i] == 0) {
	cerr << "Error: specify concentration of data " << files[i] << endl;
	exit(0);
      }
    }
    
    cerr << "==========================================" << endl;
    cerr << " Press Enter key to start" << endl;
    waitEnterKey();
    
}
  
void addConstToArray(double C, double Array[], int length) {
  for (int i=0; i<length; i++) {
    Array[i] += C;
  }
}

void copyArray(double target[], double destination[], int length) {
  for (int i=0; i<length; i++) {
    destination[i] = target[i];
  }
}

/* for debugging */
void printArr(double Y[]) {
  for (int i=0; i<lengthT; i++) {
    cout << celsius[i] << "\t" << Y[i] << endl;
  }
  exit(0);
}

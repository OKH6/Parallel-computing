// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./assignment-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>



#include <cmath>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;
double C=0;
bool CheckForColls=false;

int NumberOfBodies = 0;
int NumberOfBodiesAtStart=0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * If you need additional helper data structures, you can
 * initialise them here. Alternatively, you can introduce a
 * totally new function to initialise additional data fields and
 * call this new function from main after setUp(). Either way is
 * fine.
 *
 * This operation's semantics is not to be changed in the assignment.
 */

void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;
  C=0.01/(NumberOfBodies);
  NumberOfBodiesAtStart=NumberOfBodies;
  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


void collisions(){
  //std::cout<<"PA "<<distances[0][1]<<"PA "<<distances[1][0]<<std::endl;
  bool collis=true;
  double NMass=0;
  double pos2=0;
  while (collis==true)
  {
    collis=false;
    int pA=0,pB=0;

    for(int i=0;i<NumberOfBodies-1;i++){
      #pragma omp simd  private(pos2)
      for(int j=i+1;j<NumberOfBodies;j++){
        const double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
        );
        pos2=C*(mass[i]+mass[j]);
        if(distance<=((0.01*(mass[i]+mass[j]))/NumberOfBodiesAtStart)){
        //if(distance<=(C*(mass[i]+mass[j]))){
          //std::cout<<pos2+1e-10<<std::endl;
          pA=i;
          pB=j;
          collis=true;
        }
      }
    }
    if(collis==true){
        NMass=mass[pA]+mass[pB];
        #pragma omp simd
        for (int i = 0; i < 3; i++)
        {
          x[pA][i] = (x[pA][i] * mass[pA] + x[pB][i] * mass[pB])/NMass;
          v[pA][i] = (v[pA][i] * mass[pA] + v[pB][i] * mass[pB])/NMass;
        }
        mass[pA]=NMass;
        x[pB]=x[NumberOfBodies-1];
        v[pB]=v[NumberOfBodies-1];
        mass[pB]=mass[NumberOfBodies-1];
        NumberOfBodies--;
        //std::cout<<NumberOfBodies<<std::endl;

        

    }

    
  }
  //std::cout<<"called  "<<NumberOfBodies<<std::endl;

}




/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */
void updateBody() {
  
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();
    // FX = force along x direction
    // FY = force along y direction
    // FZ = force along z direction
  double* FX = new double[NumberOfBodies]();
  double* FY = new double[NumberOfBodies]();
  double* FZ = new double[NumberOfBodies]();
    

    //CalcDistances(distances);
  double distance=0;
  double mV=0.0;

  double* Res = new double[NumberOfBodies]();

  for(int j=0;j<NumberOfBodies;j++){
    #pragma omp simd private(distance)
    for (int i=j+1; i<NumberOfBodies; i++) {
      const double Xdiff=x[j][0]-x[i][0];
      const double Ydiff=x[j][1]-x[i][1];
      const double Zdiff=x[j][2]-x[i][2];
      distance=sqrt(Xdiff*Xdiff + Ydiff*Ydiff + Zdiff*Zdiff);
      const double Tom=(mass[i]*mass[j])/(distance *distance*distance );
      Res[i]=Tom;

    }
    for (int i = j+1; i < NumberOfBodies; i++)
    {
      const double Xdiff=x[j][0]-x[i][0];
      const double Ydiff=x[j][1]-x[i][1];
      const double Zdiff=x[j][2]-x[i][2];
      const double Fx=Xdiff*Res[i];
      const double Fy=Ydiff*Res[i];
      const double Fz=Zdiff*Res[i];
      FX[j]-=Fx;
      FY[j]-=Fy;
      FZ[j]-=Fz;
      FX[i]+=Fx;
      FY[i]+=Fy;
      FZ[i]+=Fz;
    }
  }
  
  double** fx;
  double** fv;
  fx    = new double*[NumberOfBodies];
  fv    = new double*[NumberOfBodies];

  for(int i=0;i<NumberOfBodies;i++){
    fx[i] = new double[3];
    fv[i] = new double[3];

    fx[i][0] =x[i][0] +(timeStepSize/2) * v[i][0];
    fx[i][1] =x[i][1] +(timeStepSize/2) * v[i][1];
    fx[i][2] =x[i][2] +(timeStepSize/2) * v[i][2];
    fv[i][0] =v[i][0] + (timeStepSize/2) * FX[i] / mass[i];
    fv[i][1] =v[i][1] + (timeStepSize/2) * FY[i] / mass[i];
    fv[i][2] =v[i][2] + (timeStepSize/2) * FZ[i] / mass[i];
  }
  
  for (int i = 0; i < NumberOfBodies; i++)
  {
    FX[i]=0;
    FY[i]=0;
    FZ[i]=0;
  }


  for(int j=0;j<NumberOfBodies;j++){
    #pragma omp simd reduction(min:minDx) private(distance)
    for (int i=j+1; i<NumberOfBodies; i++) {
      const double Xdiff=fx[j][0]-fx[i][0];
      const double Ydiff=fx[j][1]-fx[i][1];
      const double Zdiff=fx[j][2]-fx[i][2];
      distance=sqrt(Xdiff*Xdiff + Ydiff*Ydiff + Zdiff*Zdiff);
      const double Tom=(mass[i]*mass[j])/(distance *distance*distance );
      Res[i]=Tom;
      minDx = std::min(minDx,distance);

    }
    for (int i = j+1; i < NumberOfBodies; i++)
    {
      const double Xdiff=fx[j][0]-fx[i][0];
      const double Ydiff=fx[j][1]-fx[i][1];
      const double Zdiff=fx[j][2]-fx[i][2];
      const double Fx=Xdiff*Res[i];
      const double Fy=Ydiff*Res[i];
      const double Fz=Zdiff*Res[i];
      FX[j]-=Fx;
      FY[j]-=Fy;
      FZ[j]-=Fz;
      FX[i]+=Fx;
      FY[i]+=Fy;
      FZ[i]+=Fz;
    }
  }


  delete[] Res;
  #pragma omp simd reduction(max:maxV)
  for (int j = 0; j < NumberOfBodies; j++)
  {

    x[j][0] +=timeStepSize * fv[j][0];
    x[j][1] +=timeStepSize * fv[j][1];
    x[j][2] +=timeStepSize * fv[j][2];
    v[j][0] =v[j][0]+FX[j] / mass[j] * timeStepSize;
    v[j][1] =v[j][1]+FY[j] / mass[j] * timeStepSize;
    v[j][2] =v[j][2]+FZ[j] / mass[j] * timeStepSize;
    const double mV=std::sqrt( v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2] );
    maxV = std::max( mV,maxV );
  }
  

  for (int i = 0; i < NumberOfBodies; i++)
  {
    delete[] fx[i];
    delete[] fv[i];
  }

  delete[] fx;
  delete[] fv;
    //updateParticaleParameters(FX,FY,FZ);
  collisions();
  t += timeStepSize;
  delete[] FX;
  delete[] FY;
  delete[] FZ;
  
}



/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  closeParaviewVideoFile();

  return 0;
}
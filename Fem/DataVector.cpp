#include "DataVector.hpp"
using namespace std;

/* *** constructors *******************************************************************/
DataVector::DataVector(const Coupling &_coupling, const Vector &_values, const vectorType _type): 
                       values(_values),type(_type), coupling(_coupling)
{

}

DataVector::DataVector(const Coupling &_coupling, const int size, const vectorType _type):  
                      values(size), type(_type),coupling(_coupling)
{
    
}

DataVector::DataVector(const DataVector &rhs): coupling(rhs.coupling)
{
    type = rhs.type;
    values = rhs.values;
}

/* *** overloading operators *********************************************************/
DataVector& DataVector::operator=(const DataVector &rhs)
{
    values = rhs.values;
    type=rhs.type;
    return *this;
}

/* *** functions for conversion *****************************************************/
void DataVector::typeII_2_typeI()
{
    assert(type==typeII);
    // sum up values at cross points
    communicationCrossPoints();
    // sum up values at boundary nodes
    communicationBoundaryNodes();
    // set type
    type=typeI;
}

void DataVector::typeI_2_typeII()
{
    assert(type==typeI);

    /* *** divide values at cross points by number of processes */
  for( int k=0; k<coupling.local2globalCrossPoints.length(); k++) {
      values(k) /= coupling.crossPointsNumProcs(k+1);
  }
  /* *** divide values at boundary nodes by 2 (2=number of processes...), don't divide for cross points! */
  for (int i=1; i<=coupling.numCoupling; i++ ) {
      for (int j=2; j<=coupling.boundaryNodes[i].length()-1; j++) {            
                values(coupling.boundaryNodes[i](j)-1) /=2.;
      }
  }
    /* *** change type */
  type=typeII;
}

/* *** functions for communication ************************************************************/
void DataVector::communicationCrossPoints()
{
    
    /* *** MPI Communication for global cross points */
   Vector u_crossPoints(coupling.numCrossPoints); 
    /* *** local values at all global cross points */
   for (int i=0; i<coupling.local2globalCrossPoints.length(); i++) {
       u_crossPoints( coupling.local2globalCrossPoints(i+1)-1) = values(i);
   }
     Vector u_crossPoints_gl(coupling.numCrossPoints);
   MPI::COMM_WORLD.Allreduce(u_crossPoints.data(), u_crossPoints_gl.data(),
                             coupling.numCrossPoints, MPI::DOUBLE,MPI::SUM);

   for (int i=0 ; i<coupling.local2globalCrossPoints.length(); i++) {
       values(i) = u_crossPoints_gl(coupling.local2globalCrossPoints(i+1)-1);
   }
    
}

void DataVector::communicationBoundaryNodes()
{
    
    /* *** MPI Communication for coupling boundaries (no cross points!) */    
  for (int i=1; i<=coupling.maxColor; i++) {
       for (int j=1; j<=coupling.numCoupling; j++) {
           if (coupling.colors(j) == i && coupling.boundaryNodes[j].length()-2 >0) {
                // only communicate if there is a boundary node on coupling boundary (no cross Points!)
                int sendLength = coupling.boundaryNodes[j].length()-2;
                IndexVector sendIndex(sendLength);
                sendIndex.set(0,sendLength, coupling.boundaryNodes[j].data()+1);
                Vector u_send(sendLength),u_recv(sendLength);
                                // set local values
                for( int k=0; k<sendLength; k++) {
                    u_send(k) = values(sendIndex(k)-1);
                }
                                // get values from other processes
                MPI::COMM_WORLD.Sendrecv(u_send.data() , sendLength , MPI::DOUBLE,
                                         coupling.neighbourProcs(j)-1, 0,
                                         u_recv.data() , sendLength , MPI::DOUBLE,
                                         coupling.neighbourProcs(j)-1, 0);
                                 // add values from other processes (!! numbering is opposite !!)
                 for( int k=0; k<sendLength; k++) {
                      values( sendIndex(k)-1 ) += u_recv(sendLength-k-1);
                  }
            }
        }
  }    
}

void DataVector::resize(int N)
{
    values.resize(N);
}

/* *** producing output *************************************************************************/
void DataVector::writeData(int proc,string filename)
{
    string strproc;
  
    if(proc==0)  strproc="";
    else {
      stringstream ss;
      ss << proc;
      strproc=ss.str();
    }
    values.write(filename+strproc+".dat");
}

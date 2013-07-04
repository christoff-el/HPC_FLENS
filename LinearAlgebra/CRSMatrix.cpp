#include "CRSMatrix.hpp"

/***************************************************************************************************/
/**********************************      CRSMatrix      ********************************************/
/***************************************************************************************************/

/* *** constructors and destructor *****************************************************************/
CRSMatrix::CRSMatrix():  _rowptr(0), _colindex(0),_data(0){
}

CRSMatrix::CRSMatrix(const CRSMatrix &rhs){
		_data      = rhs.data();
		_rowptr    = rhs.rowPtr();
		_colindex  = rhs.colIndex();
		_cols      = rhs.numCols();
		_rows      = rhs.numRows();
		_nonZeros  = rhs.nonZeros();	
}

CRSMatrix::CRSMatrix(IndexVector &row, IndexVector &col, Vector &values){

	int idx, rowidx, datalen, rowdiff, offset, length=row.length();
	/* *** determine number of rows and cols */
	int nRow=0,nCol=0;
	for(int k=0;k<length;k++){
		if(row(k)>nRow) nRow=row(k);
		if(col(k)>nCol) nCol=col(k);
	}
	nRow+=1;
	nCol+=1; /* row and col index start with 0... */
	
	/* *** initialize temporary variables */
	double* data = new double[length];
	int* colindex = new int[length];
	int* rowptr = new int[nRow+1];
	int* rowoff = new int[nRow];
	memset(rowptr,0,(nRow+1)*sizeof(int));
	
	/* *** create rowptr */
	for(int k=0 ; k<length ; k++){
		rowptr[row(k)+1]++;
	}
	
	for(int k=1 ; k<nRow ; k++){
		rowptr[k+1] = rowptr[k+1]+rowptr[k];
	}
	memcpy(rowoff,rowptr,nRow*sizeof(int));
	/* *** copy column index and values */
	for(int k=0;k<length;k++){
		idx = rowoff[ row(k) ];
		rowoff[row(k) ]++;
		data[idx] = values(k);
		colindex[idx] = col(k);
	}
	/* *** sort values of each row */
	for(int k=0;k<nRow;k++){
		if(rowptr[k+1]>rowptr[k] ){
			insort(colindex+rowptr[k],data+rowptr[k],rowptr[k+1]-rowptr[k]);
		}
	}
	
	/* *** sum up entries with same indices and compress the data structures */
	idx=0; rowidx=0; rowdiff=0;
	for(int k=0;k<length;){
		while(rowptr[rowidx+1]==k){
			rowidx++;
			rowptr[rowidx]-=rowdiff;
		}
		/* *** sum upp values with similar row and column indices  */
		offset=1;
		while(k+offset < rowptr[rowidx+1] && colindex[k] == colindex[k+offset]){
			data[k] += data[k+offset];
			offset+=1;
		}
		if(data[k] != 0){
			data[idx] = data[k] ;
			colindex[idx] = colindex[k];
			idx++;
			rowdiff+=(offset-1);
		}else{
			/* *** ignore values that are 0 */
			rowdiff+=offset;
		}
		k+=offset;
	}
	datalen=idx;
	for(int k=rowidx+1;k<nRow+1;k++)
		rowptr[k]= datalen;
	
	/* *** set private variables */	
	_data = new double[datalen];
	memcpy(_data, data, datalen*sizeof(double) );
	_colindex = new int[datalen];
	memcpy(_colindex, colindex, datalen*sizeof(int) );
	_rowptr    = rowptr;
	_cols      = nCol;
	_rows      = nRow;
	_nonZeros  = datalen;
	
	delete []data;
	delete []colindex; 

}

CRSMatrix::~CRSMatrix(){
	
}

/* *** getter & setter methods ***********************************************************************/
int CRSMatrix::nonZeros() const{
	return _nonZeros;
}

int CRSMatrix::numCols() const{
	return _cols;
}

int CRSMatrix::numRows() const{
	return _rows;
}

int* CRSMatrix::rowPtr() const{
	return _rowptr;
}

int* CRSMatrix::colIndex() const{
	return _colindex;
}

double* CRSMatrix::data() const{
	return _data;
}

/* *** operators ****************************************************************************************/
CRSMatrix & CRSMatrix::operator=(const CRSMatrix &rhs){
	_data     = rhs.data();
	_rowptr   = rhs.rowPtr();
	_colindex = rhs.colIndex();
	_cols     = rhs.numCols();
	_rows     = rhs.numRows();
	_nonZeros = rhs.nonZeros();
	return *this;
}

double  CRSMatrix::operator()(const int m,const int n){
	if( (_rowptr[m+1]-_rowptr[m])>0) {
		for(int j=_rowptr[m]; j<_rowptr[m+1];j++){
			if(_colindex[j]==n) return _data[j];
		}
	}
	return 0.0;
}

double &  CRSMatrix::operator()(const int m){
	return _data[m];
}

Vector& CRSMatrix::matVec(Vector &res, Vector &x){
		assert( _cols==x.length() && _rows==res.length() );
		double tmp;
		double* xdata=x.data();
		for(int k=0; k<_rows;k++){
			tmp=0;
			for(int j=_rowptr[k];j<_rowptr[k+1];j++){
				tmp+= _data[j]*xdata[_colindex[j]];
			}
			res(k) = tmp;
		}	
		return res;
}



/* *** convert routines *********************************************************************************/
void CRSMatrix::CRS2full(Matrix & A){	
	/* *** initialize and fill  matrix */	
	A.resize(_rows,_cols);
	for (int j=0;j<_rows;j++)
		for(int k=_rowptr[j] ;k<_rowptr[j+1];k++)
			A(j,_colindex[k])=_data[k];
}

void CRSMatrix::full2CRS(Matrix &A){
	int nz = 0;     //count non-zero elements
	int idx = 0;    //data index

	/* *** count number of non-zero entries */
	for(int i=0;i<A.numRows();i++)
    for(int j=0;j<A.numCols();j++)
        if(A(i,j)!=0)
            nz++;
	/* *** initialize private variables */
	_rows = A.numRows();
	_cols = A.numCols();
	_nonZeros = nz;
	delete[] _data;
	delete []_colindex;
	delete[]_rowptr;
	_data     = new double[nz];
	_colindex = new int[nz];
	_rowptr   = new int[_rows+1];

	/* *** fill _rowptr,_colindex,_data */
	for(int i=0;i<_rows;i++){
    _rowptr[i] = idx;
    for(int j=0; j<_cols;j++){
        if(A(i,j)!=0){
            _data[idx] =  A(i,j);
            _colindex[idx]= j;
            idx++;
        }
    }
	}
	_rowptr[_rows]= nz;
}

/* *** other routines */
void CRSMatrix::writeCRS(std::string filename){
	std::fstream f;
	f.open(filename.c_str(), std::ios::out);
	f << "data:     ";
	for(int k=0;k < _nonZeros;k++) f <<_data[k]<<" ";
	f << std::endl;
	f << "colindex: ";
	for(int k=0;k<_nonZeros;k++) f <<_colindex[k]<<" ";
	f << std::endl;
	f << "rowptr:   ";
	for(int k=0;k<_rows+1;k++) f <<_rowptr[k]<<" ";
	f << std::endl;
	f.close();
}

void CRSMatrix::writeFull(std::string filename){
	/* *** convert to full and print */
	Matrix A;
	CRS2full(A);
	A.write(filename);
}

bool CRSMatrix::is_index(int row, int col){
	
	for (int j=0;j<_rows;j++)
		if(_rowptr[j]==row)
			for(int k=_rowptr[j];k<_rowptr[j+1];k++)
				if(_colindex[k]==col) return true;
	
	return false;
}

/* *** private routines *******************************************************************************************/
/* sort index and data according to the elements in index */
void CRSMatrix::insort(int * index, double* data,int length){
	int val_i; 
	double val_d;
	int j;

	for(int i=1; i<length; ++i){
		val_i = index[i]; val_d = data[i];
		j = i - 1;
		for(;;) {
			if(index[j] > val_i) {
				index[j+1] = index[j];
				data[j+1] = data[j];
				j--;
				if(j<0)
					break;
			}
			else
				break;
		}
		index[j+1] = val_i;
		data[j+1] = val_d;
	}	
}


/*******************************************************************************************************************/
/***************************************      CRSIndexMatrix      **************************************************/
/*******************************************************************************************************************/

/* *** constructors and destructor *********************************************************************************/
CRSIndexMatrix::CRSIndexMatrix(): _rowptr(0), _colindex(0),_data(0) {
}

CRSIndexMatrix::CRSIndexMatrix(const CRSIndexMatrix &rhs){
	_data     = rhs.data();
	_rowptr   = rhs.rowPtr();
	_colindex = rhs.colIndex();
	_cols     = rhs.numCols();
	_rows     = rhs.numRows();
	_nonZeros = rhs.nonZeros();
	
}


CRSIndexMatrix::CRSIndexMatrix(IndexVector &row, IndexVector &col, IndexVector &values){

		int idx, rowidx, datalen, rowdiff, offset, length=row.length();
		/* *** determine number of rows and cols */
		int nRow=0,nCol=0;
		for(int k=0;k<length;k++){
			if(row(k)>nRow) nRow=row(k);
			if(col(k)>nCol) nCol=col(k);
		}
		nRow+=1;
		nCol+=1; /* row and col index start with 0... */

		/* *** initialize temporary variables */
		int* data = new int[length];
		int* colindex = new int[length];
		int* rowptr = new int[nRow+1];
		int* rowoff = new int[nRow];
		memset(rowptr,0,(nRow+1)*sizeof(int));

		/* *** create rowptr */
		for(int k=0 ; k<length ; k++){
			rowptr[row(k)+1]++;
		}

		for(int k=1 ; k<nRow ; k++){
			rowptr[k+1] = rowptr[k+1]+rowptr[k];
		}
		memcpy(rowoff,rowptr,nRow*sizeof(int));
		/* *** copy column index and values */
		for(int k=0;k<length;k++){
			idx = rowoff[ row(k) ];
			rowoff[row(k) ]++;
			data[idx] = values(k);
			colindex[idx] = col(k);
		}
		/* *** sort values of each row */
		for(int k=0;k<nRow;k++){
			if(rowptr[k+1]>rowptr[k] ){
				insort(colindex+rowptr[k],data+rowptr[k],rowptr[k+1]-rowptr[k]);
			}
		}

		/* *** sum up entries with same indices and compress the data structures */
		idx=0; rowidx=0; rowdiff=0;
		for(int k=0;k<length;){
			while(rowptr[rowidx+1]==k){
				rowidx++;
				rowptr[rowidx]-=rowdiff;
			}
			/* *** sum upp values with similar row and column indices  */
			offset=1;
			while(k+offset < rowptr[rowidx+1] && colindex[k] == colindex[k+offset]){
				data[k] += data[k+offset];
				offset+=1;
			}
			if(data[k] != 0){
				data[idx] = data[k] ;
				colindex[idx] = colindex[k];
				idx++;
				rowdiff+=(offset-1);
			}else{
				/* *** ignore values that are 0 */
				rowdiff+=offset;
			}
			k+=offset;
		}
		datalen=idx;
		for(int k=rowidx+1;k<nRow+1;k++)
			rowptr[k]= datalen;

		/* *** set private variables */	
		_data = new int[datalen];
		memcpy(_data, data, datalen*sizeof(int) );
		_colindex = new int[datalen];
		memcpy(_colindex, colindex, datalen*sizeof(int) );
		_rowptr    = rowptr;
		_cols      = nCol;
		_rows      = nRow;
		_nonZeros  = datalen;

		delete []data;
		delete []colindex; 

	}

CRSIndexMatrix::~CRSIndexMatrix(){
}

/* *** getter & setter methods *****************************************************************************/
// int CRSIndexMatrix::get(int row,int col) const{
// 	
// 	if( (_rowptr[row+1]-_rowptr[row])>0) {
// 		for(int j=_rowptr[row]; j<_rowptr[row+1];j++){
// 			if(_colindex[j]==col) return _data[j];
// 		}
// 	}
// 	return 0;
// }

// int CRSIndexMatrix::get(int k) const{
// 	return _data[k];
// }

int CRSIndexMatrix::nonZeros() const{
	return _nonZeros;
}

int CRSIndexMatrix::numCols() const{
	return _cols;
}

int CRSIndexMatrix::numRows() const{
	return _rows;
}

int* CRSIndexMatrix::rowPtr() const{
	return _rowptr;
}

int* CRSIndexMatrix::colIndex() const{
	return _colindex;
}

int* CRSIndexMatrix::data() const{
	return _data;
}

/* *** operators *******************************************************************************************/
CRSIndexMatrix& CRSIndexMatrix::operator=(const CRSIndexMatrix &rhs){
	_data=rhs.data();
	_rowptr=rhs.rowPtr();
	_colindex=rhs.colIndex();
	_cols = rhs.numCols();
	return *this;
}
int  CRSIndexMatrix::operator()(const int m,const int n){
	if( (_rowptr[m+1]-_rowptr[m])>0) {
		for(int j=_rowptr[m]; j<_rowptr[m+1];j++){
			if(_colindex[j]==n) return _data[j];
		}
	}
	return 0;
}

int&  CRSIndexMatrix::operator()(const int m){
	return _data[m];
}

/* *** convert routines ***********************************************************************************/
void CRSIndexMatrix::CRS2full(IndexMatrix& A){	
	/* *** initialize and fill  matrix */	
	A.resize(_rows,_cols);
	for (int j=0;j<_rows;j++)
		for(int k=_rowptr[j] ;k<_rowptr[j+1];k++)
			A(j,_colindex[k])=_data[k];
}

void CRSIndexMatrix::full2CRS(IndexMatrix &A){
	int nz = 0;     //count non-zero elements
	int idx = 0;    //data index

	/* *** count number of non-zero entries */
	for(int i=0;i<A.numRows();i++)
    for(int j=0;j<A.numCols();j++)
        if(A(i,j)!=0)
            nz++;
	/* *** initialize private variables */
	_rows = A.numRows();
	_cols = A.numCols();
	_nonZeros = nz;
	delete[] _data;
	delete []_colindex;
	delete[]_rowptr;
	_data     = new int[nz];
	_colindex = new int[nz];
	_rowptr   = new int[_rows+1];

	/* *** fill _rowptr,_colindex,_data */
	for(int i=0;i<_rows;i++){
    _rowptr[i] = idx;
    for(int j=0; j<_cols;j++){
        if(A(i,j)!=0){
            _data[idx] =  A(i,j);
            _colindex[idx]= j;
            idx++;
        }
    }
	}
	_rowptr[_rows]= nz;
}

/* *** other routines *******************************************************************************************/
void CRSIndexMatrix::writeCRS(std::string filename){
	std::fstream f;
	f.open(filename.c_str(), std::ios::out);
	f << "data:     ";
	for(int k=0;k < _nonZeros;k++) f <<_data[k]<<" ";
	f << std::endl;
	f << "colindex: ";
	for(int k=0;k<_nonZeros;k++) f <<_colindex[k]<<" ";
	f << std::endl;
	f << "rowptr:   ";
	for(int k=0;k<_rows+1;k++) f <<_rowptr[k]<<" ";
	f << std::endl;
	f.close();
}

void CRSIndexMatrix::writeFull(std::string filename){
	/* *** convert to full and print */
	IndexMatrix A;
	CRS2full(A);
	A.write(filename);
}

bool CRSIndexMatrix::is_index(int row, int col){
	
	for (int j=0;j<_rows;j++)
		if(_rowptr[j]==row)
			for(int k=_rowptr[j];k<_rowptr[j+1];k++)
				if(_colindex[k]==col) return true;
	
	return false;
}

/* *** private routines ******************************************************************************************/
/* sort index and data according to the elements in index */
void CRSIndexMatrix::insort(int * index, int* data,int length){
	int val_i; 
	int val_d;
	int j;

	for(int i=1; i<length; ++i){
		val_i = index[i]; val_d = data[i];
		j = i - 1;
		for(;;) {
			if(index[j] > val_i) {
				index[j+1] = index[j];
				data[j+1] = data[j];
				j--;
				if(j<0)
					break;
			}
			else
				break;
		}
		index[j+1] = val_i;
		data[j+1] = val_d;
	}	
}

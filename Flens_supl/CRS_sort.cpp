#ifndef CRS_SORT_CPP
#define CRS_SORT_CPP 1

#include "CRS_sort.h"

namespace flens{

void
CRS_sort(IVector &_data, IVector &row, IVector &col, IVector &values)
{
	//--- Typedefs for FLENS
    // Used to set values from a C-Array
    typedef flens::DenseVector<flens::Array<int>::View> 	VecView;
    typedef flens::Array<int>::View 						View;
    //--- Typedefs for FLENS


	int idx, rowidx, rowdiff, offset, length=row.length();
	/* *** determine number of rows and cols */
	int nRow=0,nCol=0;
	for(int k=1; k<=length; ++k)
	{
		if(row(k)>nRow) nRow=row(k);
		if(col(k)>nCol) nCol=col(k);
	}
	
	/* *** initialize temporary variables */
	int* data 		= new int[length];
	int* colindex 	= new int[length];
	int* rowptr 	= new int[nRow+1];
	int* rowoff 	= new int[nRow];
	memset(rowptr,0,(nRow+1)*sizeof(int));

	/* *** create rowptr */
	for(int k=1 ; k<=length; ++k)
	{
		++rowptr[row(k)];
	}
	
	for(int k=1 ; k<nRow ; ++k)
	{
		rowptr[k+1] = rowptr[k+1]+rowptr[k];
	}

	memcpy(rowoff,rowptr,nRow*sizeof(int));
	/* *** copy column index and values */
	for(int k=1; k<=length; ++k)
	{
		idx = rowoff[row(k)-1];
		++rowoff[row(k)-1];
		data[idx] = values(k);
		colindex[idx] = col(k);
	}
	/* *** sort values of each row */
	for(int k=0; k<nRow; ++k)
	{
		if(rowptr[k+1]>rowptr[k])
		{
			CRS_insort(colindex+rowptr[k], data+rowptr[k],rowptr[k+1]-rowptr[k]);
		}
	}
	
	/* *** sum up entries with same indices and compress the data structures */
	idx=0; rowidx=0; rowdiff=0;
	for(int k=0; k<length;)
	{
		while(rowptr[rowidx+1]==k)
		{
			rowidx++;
			rowptr[rowidx]-=rowdiff;
		}
		/* *** sum upp values with similar row and column indices  */
		offset=1;
		while(k+offset<rowptr[rowidx+1] && colindex[k] == colindex[k+offset])
		{
			data[k] += data[k+offset];
			offset+=1;
		}
		if(data[k] != 0)
		{
			data[idx] = data[k] ;
			colindex[idx] = colindex[k];
			++idx;
			rowdiff += (offset-1);
		}
		else
		{
			/* *** ignore values that are 0 */
			rowdiff+=offset;
		}
		k+=offset;
	}
	
	/* *** set private variables */
	_data.resize(idx);
	_data = VecView(View(idx, data));
	
	delete[] data;
	delete[] colindex;
	delete[] rowptr;
	delete[] rowoff;
};

void
CRS_insort(int * index, int* data, int length)
{
	int val_i; 
	int val_d;
	int j;

	for(int i=1; i<length; ++i)
	{
		val_i = index[i]; val_d = data[i];
		j = i - 1;
		for(;;)
		{
			if(index[j] > val_i)
			{
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
};

};// namespace flens

#endif //CRS_SORT_CPP
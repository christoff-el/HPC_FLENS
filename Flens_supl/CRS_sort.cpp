#ifndef CRS_SORT_CPP
#define CRS_SORT_CPP 1

#include "CRS_sort.h"

namespace flens{

template <typename DVector, typename IVector>
void
CRS_sort(DVector &_data, IVector &row, IVector &col, DVector &values)
{
	//--- Typedefs for FLENS
	typedef typename DVector::ElementType									ElementType;
	typedef typename IVector::ElementType									IndexType;
    typedef flens::DenseVector<typename flens::Array<ElementType>::View> 	VecView;
    typedef typename flens::Array<ElementType>::View 						View;
    //--- Typedefs for FLENS


	IndexType idx, rowidx, rowdiff, offset, length=row.length();
	/* *** determine number of rows and cols */
	IndexType nRow=0,nCol=0;
	for(IndexType k=1; k<=length; ++k)
	{
		if(row(k)>nRow) nRow=row(k);
		if(col(k)>nCol) nCol=col(k);
	}
	
	/* *** initialize temporary variables */
	ElementType* data 		= new ElementType[length];
	IndexType* colindex 	= new IndexType[length];
	IndexType* rowptr 		= new IndexType[nRow+1];
	IndexType* rowoff 		= new IndexType[nRow];
	memset(rowptr,0,(nRow+1)*sizeof(IndexType));

	/* *** create rowptr */
	for(IndexType k=1 ; k<=length; ++k)
	{
		++rowptr[row(k)];
	}
	
	for(IndexType k=1 ; k<nRow ; ++k)
	{
		rowptr[k+1] = rowptr[k+1]+rowptr[k];
	}

	memcpy(rowoff,rowptr,nRow*sizeof(IndexType));
	/* *** copy column index and values */
	for(IndexType k=1; k<=length; ++k)
	{
		idx = rowoff[row(k)-1];
		++rowoff[row(k)-1];
		data[idx] = values(k);
		colindex[idx] = col(k);
	}
	/* *** sort values of each row */
	for(IndexType k=0; k<nRow; ++k)
	{
		if(rowptr[k+1]>rowptr[k])
		{
			CRS_insort(colindex+rowptr[k], data+rowptr[k],rowptr[k+1]-rowptr[k]);
		}
	}
	
	/* *** sum up entries with same indices and compress the data structures */
	idx=0; rowidx=0; rowdiff=0;
	for(IndexType k=0; k<length;)
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

template <typename Data, typename Index>
void
CRS_insort(Index* index, Data* data, Index length)
{
	Index val_i; 
	Data val_d;
	Index j;

	for(Index i=1; i<length; ++i)
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

#endif	//CRS_SORT_CPP

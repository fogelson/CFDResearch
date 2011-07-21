/*
 * MacroscaleObjects.h
 *
 *  Created on: Jun 19, 2011
 *      Author: fogelson
 */

#ifndef MACROSCALEOBJECTS_H_
#define MACROSCALEOBJECTS_H_

#include <blitz/array.h>

using namespace blitz;

namespace CFD{
	template<int N_rank>
	class SymmetricTensorArrayOrder : public GeneralArrayStorage<N_rank>{
	private:
	    typedef GeneralArrayStorage<N_rank> T_base;
	    typedef _bz_typename T_base::noInitializeFlag noInitializeFlag;
	    using T_base::ordering_;
	    using T_base::ascendingFlag_;
	    using T_base::base_;
	public:
	    SymmetricTensorArrayOrder() : GeneralArrayStorage<N_rank>(noInitializeFlag()){
	    	if(N_rank != 3){
	    		std::cout << "Error. Any Blitz array with the SymmetricTensorArrayOrder storage order"
	    				<< " should be three-dimensional." << std::endl;
	    	}
	        for(int i = 0; i < N_rank; i++){
	        	ordering_(i) = N_rank - i - 1;
	        }
	    	ascendingFlag_ = true;
	        base_ = 0;
	    }
	};

	template<int N_rank>
	class VelocityArrayOrder : public GeneralArrayStorage<N_rank>{
	private:
		typedef GeneralArrayStorage<N_rank> T_base;
	    typedef _bz_typename T_base::noInitializeFlag noInitializeFlag;
	    using T_base::ordering_;
	    using T_base::ascendingFlag_;
	    using T_base::base_;
	public:
	    VelocityArrayOrder() : GeneralArrayStorage<N_rank>(noInitializeFlag()){
	    	if(N_rank != 3){
	    		std::cout << "Error. Any Blitz array with the VelocityArrayOrder storage order"
	    				<< " should be three-dimensional." << std::endl;
	    	}
	    	for(int i = 0; i < N_rank; i++){
	    		ordering_(i) = N_rank - i - 1;
	    	}
	    	ascendingFlag_ = true;
	    	base_ = 0;
	    }
	};

	class VelocityArray : Array<double,3>{
		int N;
	public:
		VelocityArray(double *U, int N){
			VelocityArrayOrder<3> storage;
			setStorage(storage);
			this->N = N;
		}
	};
}


#endif /* MACROSCALEOBJECTS_H_ */

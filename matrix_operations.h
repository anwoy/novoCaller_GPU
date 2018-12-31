#ifndef GAURD_matrix_operations
#define GAURD_matrix_operations

namespace std {
// determinant and inverse
// *************

int determinant_sign(const boost::numeric::ublas::permutation_matrix<std ::size_t> &pm){
	int pm_sign=1;
	size_t size = pm.size();
	for (size_t i = 0; i < size; ++i)
		if (i != pm(i))
			pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
	return pm_sign;
}
 
double determinant( boost::numeric::ublas::matrix<double> &m ){
	boost::numeric::ublas::permutation_matrix<size_t> pm(m.size1());
	double det = 1.0;
	if( boost::numeric::ublas::lu_factorize(m,pm) ) {
		det = 0.0;
	}
	else {
		for(int i = 0; i < m.size1(); i++)
			det *= m(i,i); // multiply by elements on diagonal
		det = det * determinant_sign( pm );
	}
	return det;
}

double determinant(const vector<vector<double> > &cov){
	if(cov.size() != cov[0].size())
		throw(domain_error("ERROR1 in determinant"));
	const size_t N(cov.size());
	boost::numeric::ublas::matrix<double> M(N,N);
	for(size_t i(0);i<N;++i)
		for(size_t j(0);j<N;++j)
			M(i,j) = cov[i][j];
	return determinant(M);
}

template<class T>
bool InvertMatrix(const boost::numeric::ublas::matrix<T> &input, boost::numeric::ublas::matrix<T> &inverse){
	
	// create a working copy of the input
	boost::numeric::ublas::matrix<T> A(input);
	
	// create a permutation matrix for the LU-factorization
	boost::numeric::ublas::permutation_matrix<size_t> pm(A.size1());
	
	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;
	
	// create identity matrix of "inverse"
	inverse.assign(boost::numeric::ublas::identity_matrix<T> (A.size1()));
	
	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
	
	return true;
}


vector<vector<double> > InvertMatrix(const vector<vector<double> > &V){
	if(V.size() != V[0].size())
		throw(domain_error("ERROR1 in InvertMatrix"));
	const size_t N(V.size());
	boost::numeric::ublas::matrix<double> M(N,N),inverse(N,N);
	for(size_t i(0);i<N;++i)
		for(size_t j(0);j<N;++j)
			M(i,j)=V[i][j];
	InvertMatrix(M,inverse);
	
	vector<vector<double> > ans(N,vector<double>(N));
	for(size_t i(0);i<N;++i)
		for(size_t j(0);j<N;++j)
			ans[i][j] = inverse(i,j);
	return ans;
}

// *************

// extraction
// *************


vector<double> extract(const vector<double> &V,const vector<bool> &index){
	if(index.size()!=V.size())
		throw(domain_error("ERROR1 in extract"));
	size_t count(0);
	for(size_t i(0);i<index.size();++i)
		if(index[i])
			count += 1;
	vector<double> ans(count);
	count=0;
	for(size_t i(0);i<index.size();++i){
		if(index[i]){
			ans[count]=V[i];
			++count;
		}
	}
	return ans;
}


vector<vector<double> > extract(const vector<vector<double> > &V,const vector<bool> &index1,const vector<bool> &index2){
	if(V.size()!=index1.size() || V[0].size()!=index2.size())
		throw(domain_error("ERROR2 in extract"));
	size_t count1(0),count2(0);
	for(size_t i(0);i<index1.size();++i)
		if(index1[i])
			count1 += 1;
	for(size_t i(0);i<index2.size();++i)
		if(index2[i])
			count2 += 1;
		
	vector<vector<double> > ans(count1,vector<double>(count2));
	count1=0;
	for(size_t i(0);i<index1.size();++i){
		if(index1[i]){
			count2=0;
			for(size_t j(0);j<index2.size();++j){
				if(index2[j]){
					ans[count1][count2] = V[i][j];
					count2 += 1;
				}
			}
			count1 += 1;
		}
	}
	return ans;
}

// *************


// log and exponent
// *************



vector<double> LOG_1D(const vector<double> &V){
	vector<double> ans(V.size());
	for(size_t i(0);i<V.size();++i)
		ans[i] = log(V[i]);
	return ans;
}

void LOG_1D_INPLACE(vector<double> &V){
	for(size_t i(0);i<V.size();++i)
		V[i] = log(V[i]);
}


vector<vector<double> > LOG_2D(const vector<vector<double> > &V){
	vector<vector<double> > ans(V.size(),vector<double>(V[0].size()));
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[0].size();++j)
			ans[i][j] = log( V[i][j] );
	return ans;
}

void LOG_2D_INPLACE(vector<vector<double> > &V){
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[i].size();++j)
			V[i][j] = log(V[i][j]);
}


vector<double> EXP_1D(const vector<double> &V){
	vector<double> ans(V.size());
	for(size_t i(0);i<V.size();++i)
		ans[i] = exp(V[i]);
	return ans;
}

void EXP_1D_INPLACE(vector<double> &V){
	for(size_t i(0);i<V.size();++i)
		V[i] = exp(V[i]);
}

vector<vector<double> > EXP_2D(const vector<vector<double> > &V){
	vector<vector<double> > ans(V.size(),vector<double>(V[0].size()));
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[0].size();++j)
			ans[i][j] = exp( V[i][j] );
	return ans;
}

void EXP_2D_INPLACE(vector<vector<double> > &V){
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[i].size();++j)
			V[i][j] = exp(V[i][j]);
}

// *************


// EYE
// *************

vector<vector<double> > EYE(const size_t &D){
	vector<vector<double> > ans(D,vector<double>(D,0));
	for(size_t i(0);i<D;++i)
		ans[i][i]=1;
	return ans;
}

// *************


// matrix multiplication
// *************

vector<double> DOT(const vector<double> &V,const vector<vector<double> > &M){
	if(V.size()!=M.size())
		throw(domain_error("ERROR1 in DOT"));
	
	vector<double> ans(M[0].size());
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += V[j]*M[j][i];
	}
	return ans;
}

void DOT(const vector<double> &V,const vector<vector<double> > &M,vector<double> &ans){
	if(V.size()!=M.size())
		throw(domain_error("ERROR12 in DOT"));
	if(M[0].size() != ans.size())
		throw(domain_error("ERROR123 in DOT"));
	
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += V[j]*M[j][i];
	}
}


vector<double> DOT(const vector<vector<double> > &M,const vector<double> &V){
	if(M[0].size()!=V.size())
		throw(domain_error("ERROR2 in DOT"));
	vector<double> ans(M.size());
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += M[i][j]*V[j];
	}
	return ans;
}


void DOT(const vector<vector<double> > &M,const vector<double> &V,vector<double> &ans){
	if(M[0].size()!=V.size())
		throw(domain_error("ERROR21 in DOT"));
	if(M.size() != ans.size())
		throw(domain_error("ERROR213 in DOT"));
	for(size_t i(0);i<ans.size();++i){
		ans[i]=0;
		for(size_t j(0);j<V.size();++j)
			ans[i] += M[i][j]*V[j];
	}
}

double DOT(const vector<double> &v1,const vector<double> &v2){
	if(v1.size() != v2.size())
		throw(domain_error("ERROR3 in DOT"));
	double ans(0);
	for(size_t i(0);i<v1.size();++i)
		ans += v1[i]*v2[i];
	return ans;
}

vector<vector<double> > DOT(const vector<vector<double> > &v1,const vector<vector<double> > &v2){
	if(v1[0].size()!=v2.size())
		throw(domain_error("ERROR4 in DOT"));
	vector<vector<double> > ans(v1.size(),vector<double>(v2[0].size()));
	for(size_t i(0);i<ans.size();++i){
		for(size_t j(0);j<ans[0].size();++j){
			ans[i][j] = 0;
			for(size_t k(0);k<v1[0].size();++k)
				ans[i][j] += v1[i][k]*v2[k][j];
		}
	}
	return ans;
}

void DOT(const vector<vector<double> > &v1,const vector<vector<double> > &v2,vector<vector<double> > &ans){
	if(v1[0].size()!=v2.size())
		throw(domain_error("ERROR5 in DOT"));
	if(ans.size() != v1.size() || ans[0].size() != v2[0].size())
		throw(domain_error("ERROR6 in DOT"));
	for(size_t i(0);i<ans.size();++i){
		for(size_t j(0);j<ans[0].size();++j){
			ans[i][j] = 0;
			for(size_t k(0);k<v1[0].size();++k)
				ans[i][j] += v1[i][k]*v2[k][j];
		}
	}
}

// *************


// matrix addition
// *************


vector<double> ADD(const vector<double> &v1,const vector<double> &v2){
	if(v1.size() != v2.size())
		throw(domain_error("ERROR1 in ADD"));
	vector<double> ans(v1.size());
	for(size_t i(0);i<v1.size();++i)
		ans[i] = v1[i]+v2[i];
	return ans;
}

void ADD(const vector<double> &v1,const vector<double> &v2,vector<double> &ans){
	if(v1.size() != v2.size() || v1.size() != ans.size())
		throw(domain_error("ERROR2 in ADD"));
	for(size_t i(0);i<v1.size();++i)
		ans[i] = v1[i]+v2[i];
}

vector<vector<double> > ADD(const vector<vector<double> > &v1,const vector<vector<double> > &v2){
	if(v1.size() != v2.size() || v1[0].size() != v2[0].size())
		throw(domain_error("ERROR3 in ADD"));
	vector<vector<double> > ans(v1.size(),vector<double>(v1[0].size()));
	for(size_t i(0);i<v1.size();++i)
		for(size_t j(0);j<v1[0].size();++j)
			ans[i][j] = v1[i][j]+v2[i][j];
	return ans;
}

void ADD(const vector<vector<double> > &v1,const vector<vector<double> > &v2,vector<vector<double> > &ans){
	if(v1.size() != v2.size() || v1[0].size() != v2[0].size() || v1.size() != ans.size() || v1[0].size() != ans[0].size())
		throw(domain_error("ERROR4 in ADD"));
	for(size_t i(0);i<v1.size();++i)
		for(size_t j(0);j<v1[0].size();++j)
			ans[i][j] = v1[i][j]+v2[i][j];
}

// *************

// scalar multiplication with vectors/matrices (no need for division)
// *************

void operator*=(vector<double> &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		V[i] = S*V[i];
}

void operator*=(vector<vector<double> > &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[0].size();++j)
			V[i][j] = S*V[i][j];
}

// *************

// scalar addition with vectors/matrices (no need for subtraction)
// *************

void operator+=(vector<double> &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		V[i] = S+V[i];
}

void operator+=(vector<vector<double> > &V,const double &S){
	for(size_t i(0);i<V.size();++i)
		for(size_t j(0);j<V[0].size();++j)
			V[i][j] = S+V[i][j];
}

// *************

// sum of 1D and 2D matrices
// *************
double SUM(const vector<double> &V){
	double ans(0);
	for(size_t i(0);i<V.size();++i)
		ans += V[i];
	return ans;
}

vector<double> SUM(const vector<vector<double> > &V,const size_t &dim){
	if(dim > 1)
		throw(domain_error("ERROR1 in vector<double> SUM(const vector<vector<double> > &V,const size_t &dim)"));
	
	if(dim == 0){
		vector<double> ans(V[0].size(),0);
		for(size_t i(0);i<V.size();++i)
			for(size_t j(0);j<V[0].size();++j)
				ans[j] += V[i][j];
		return ans;
	}
	else{
		vector<double> ans(V.size(),0);
		for(size_t i(0);i<V.size();++i)
			for(size_t j(0);j<V[0].size();++j)
				ans[i] += V[i][j];
		return ans;
	}
}

// *******************


// copying of vectors and matricers
// *************

void COPY(const vector<double> &v1,vector<double> &v2){
	if(v2.size() != v1.size())
		v2 = v1;
	else
		for(size_t i(0);i<v2.size();++i)
			v2[i] = v1[i];
}

void COPY(const vector<vector<double> > &v1,vector<vector<double> > &v2){
	if(v2.size() != v1.size())
		v2 = v1;
	else
		for(size_t i(0);i<v2.size();++i)
			COPY(v1[i],v2[i]);
}

void COPY(const vector<vector<vector<double> > > &v1,vector<vector<vector<double> > > &v2){
	if(v2.size() != v1.size())
		v2 = v1;
	else
		for(size_t i(0);i<v2.size();++i)
			COPY(v1[i],v2[i]);
}



// *******************


// MAX and ARGMAX

double MAX(const vector<double> &V){
	if(V.size() > 0){
		double maxval( V[0] );
		for(size_t i(0);i<V.size();++i)
			if(V[i] > maxval)
				maxval = V[i];
		return maxval;
	}
	else
		throw(domain_error("ERROR1 in MAX"));
}


size_t ARGMAX(const vector<double> &V){
	if(V.size() > 0){
		size_t argmax(0);
		double maxval(V[0]);
		for(size_t i(0);i<V.size();++i)
			if(V[i] > maxval){
				maxval = V[i];
				argmax = i;
			}
		return argmax;
	}
	else
		throw(domain_error("ERROR1 in ARGMAX"));
}

}

#endif































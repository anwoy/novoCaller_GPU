#ifndef GAURD_defs01
#define GAURD_defs01

#define NUMBER_OF_PREFIXES 9
#define LOG_2 0.6931471805599453
#define A_prior 1.01
#define B_prior 1.01
#define D_prior 1.01
#define THRESHOLD 1e-10
#define MAX_ITER 500
#define MUTATION_RATE 1e-11

////////////////
struct initialization_info {
	int CHROM_col, POS_col, ID_col, REF_col, ALT_col, QUAL_col, FILTER_col, INFO_col, FORMAT_col;
	int sample_count;
	
	initialization_info(std::string &);
	
	void set_header_info(const std::string &);
};

std::ostream & operator << (std::ostream & out, const initialization_info & O) {
	using namespace std;
	out<< "initialization_info:\n";
	out<< "********************\n";
	cout<< "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"<<"\n";
	cout<< O.CHROM_col <<"\t"<< O.POS_col <<"\t"<< O.ID_col <<"\t"<< O.REF_col <<"\t"<< O.ALT_col <<"\t"<< O.QUAL_col <<"\t"<< O.FILTER_col <<"\t"<< O.INFO_col <<"\t"<< O.FORMAT_col <<"\n";
	cout<< "sample_count = " << O.sample_count <<"\n";
	
	return out;
}

initialization_info::initialization_info(std::string & line) {
	set_header_info(line);
}


void initialization_info::set_header_info(const std::string & header_line) {
	using namespace std;
	CHROM_col = POS_col = ID_col = REF_col = ALT_col = QUAL_col = FILTER_col = INFO_col = FORMAT_col = -1;
	const vector<string> SL( split(header_line) );
	for(size_t i(0) ; i < SL.size() ; ++i) {
		if(SL[i] == "#CHROM")
			CHROM_col  = i;
		if(SL[i] == "POS")
			POS_col = i;
		if(SL[i] == "ID")
			ID_col = i;
		if(SL[i] == "REF")
			REF_col = i;
		if(SL[i] == "ALT")
			ALT_col = i;
		if(SL[i] == "QUAL")
			QUAL_col = i;
		if(SL[i] == "FILTER")
			FILTER_col = i;
		if(SL[i] == "INFO")
			INFO_col = i;
		if(SL[i] == "FORMAT")
			FORMAT_col = i;
	}
	sample_count = SL.size() - NUMBER_OF_PREFIXES;
}

////////////////

struct AD_info {
	thrust::host_vector<int> AD0, AD1;
	thrust::device_vector<int> DAD0,DAD1;
	int FORMAT_length, AD_position, CHROM, POS;
	
	AD_info(const initialization_info &);
	
	bool read_once(std::istream &, const initialization_info &);
	
	void FORMAT_info_extract(const std::string &);
	void CHROM_extract(const std::string &);
	std::string integer_to_chromosome() const;
	void POS_extract(const std::string &);
};

std::string AD_info::integer_to_chromosome() const {
	using namespace std;
	string temp;
	if(CHROM == 0)
		temp = string("M");
		else
			if(CHROM == 23)
				temp = string("X");
			else
				if(CHROM == 24)
					temp = string("Y");
				else
					temp = myitoa(CHROM);
	
	const string ans(string("chr") + temp);
	return ans;
}

void AD_info::POS_extract(const std::string & S) {
	using namespace std;
	POS = atoi( S.c_str() );
}


void AD_info::CHROM_extract(const std::string & S) {
	using namespace std;
	string S_temp, S2;
	if(S.size() > 3)
		S_temp = string(S.begin(), S.begin() + 3);
	else
		S_temp=S;
	
	if(S_temp == "chr")
		S2 = string(S.begin() + 3, S.end());
	else
		S2=S;
	
	const int num( atoi( S2.c_str() ) );
	if(num >= 1 and num <= 22)
		CHROM  = num;
	else
		if(S2 == "M" or S2 == "m")
			CHROM = 0;
		else
			if(S2 == "X" or S2 == "x")
				CHROM = 23;
			else
				if(S2 == "Y" or S2 == "y")
					CHROM = 24;
				else
					throw(domain_error("ERROR in CHROM_extract"));
}

std::ostream & operator << (std::ostream & out, const AD_info & O) {
	using namespace std;
	out<< "AD_info:\n";
	out<< "********\n";
	cout<< "AD_position = " << O.AD_position <<"\n";
	cout<< "FORMAT_length = " << O.FORMAT_length <<"\n";
	cout<< "CHROM, POS = "<< O.CHROM <<", "<< O.POS <<"\n";
	if(O.AD0.size() != O.AD1.size() or O.DAD0.size() != O.DAD1.size() or O.AD0.size() != O.DAD0.size() )
		throw(domain_error("ERROR in operator << (std::ostream & out, const AD_info & O)"));
	for(int i(0) ; i < O.AD0.size() ; ++i)
		cout<< i << ")\t" << O.AD0[i] <<"\t"<< O.AD1[i] <<"\t"<< O.DAD0[i] <<"\t"<< O.DAD1[i] <<"\n";
	return out;
}

AD_info::AD_info(const initialization_info & IIO) {
	using namespace thrust;
	AD0 = host_vector<int>(IIO.sample_count);
	AD1 = host_vector<int>(IIO.sample_count);
	DAD0 = device_vector<int>(IIO.sample_count);
	DAD1 = device_vector<int>(IIO.sample_count);
}

void AD_info::FORMAT_info_extract(const std::string & FORMAT) {
	using namespace std;
	AD_position = FORMAT_length = -1;
	const vector<string> SL( split(FORMAT, ":") );
	FORMAT_length = SL.size();
	for(int i(0) ; i < SL.size() ; ++i)
		if(SL[i] == "AD")
			AD_position = i;
	if(AD_position == -1 or FORMAT_length == -1)
		throw(domain_error("ERROR in AD_info::AD_info_extract"));
}

bool AD_info::read_once(std::istream & fin, const initialization_info & IIO) {
	std::string temp;
	for(int i(0) ; i < NUMBER_OF_PREFIXES ; ++i) {
		fin >> temp;
		if(not fin)
			return false;
		if(i == IIO.FORMAT_col)
			FORMAT_info_extract(temp);
		if(i == IIO.CHROM_col)
			CHROM_extract(temp);
		if(i == IIO.POS_col)
			POS_extract(temp);
	}
	for(int i(0) ; i < AD0.size() ; ++i) {
		fin >> temp;
		std::vector<std::string> SL( split(temp, ":") );
		if(SL.size() != FORMAT_length) {
			AD0[i] = 0;
			AD1[i] = 0;
			continue;
		}
		temp = SL[AD_position];
		SL = split(temp, ",");
		if(SL.size() <= 1) {
			AD0[i] = 0;
			AD1[i] = 0;
			continue;
		}
		AD0[i] = atoi( SL[0].c_str() );
		AD1[i] = 0;
		for(int j(1) ; j < SL.size() ; ++j)
			AD1[i] += atoi( SL[j].c_str() );
	}
	
	thrust::copy(AD0.begin(), AD0.end(), DAD0.begin());
	thrust::copy(AD1.begin(), AD1.end(), DAD1.begin());
	
	return true;
}

////////////////

struct sample_info {
	std::vector<int> control_index;
	std::vector<std::vector<int> > trios;
	std::vector<std::string> sample_names;
	
	thrust::device_vector<int> device_control_index;
	
	thrust::device_vector<int> device_P1_index;
	thrust::device_vector<int> device_P2_index;
	thrust::device_vector<int> device_child_index;
	
	sample_info(const std::string & header_line, const std::string & control_filename, const std::string & trio_filename);
	
	int find_in_line(const std::string &, const std::vector<std::string> &);
};


std::ostream & operator << (std::ostream & out, const sample_info & O) {
	using namespace std;
	out<< "sample_info:\n";
	out<< "************\n";
	out<< "control_index:\n";
	out<< O.control_index <<"\t:"<< "number_of_control = "<< O.control_index.size() <<"\n";
	out<< "device_control_index:\n";
	out<< O.device_control_index <<"\t:"<< "number_of_control = "<< O.device_control_index.size() <<"\n";
	out<< "trios:\n";
	for(int i(0);i<O.trios.size();++i)
		out<< O.trios[i] <<"\n";
	out<< "number of trios = " << O.trios.size() <<"\n";
	out<< "device_trio_indexes:\n";
	out<< O.device_P1_index    <<"\t"<< "length = " << O.device_P1_index.size()    <<"\n";
	out<< O.device_P2_index    <<"\t"<< "length = " << O.device_P2_index.size()    <<"\n";
	out<< O.device_child_index <<"\t"<< "length = " << O.device_child_index.size() <<"\n";
	for(int i(0);i<O.sample_names.size();++i)
		out<< i <<") "<< O.sample_names[i] <<"\t";
	out<< ":no_of_samples = "<< O.sample_names.size() <<"\n";
	return out;
}


int sample_info::find_in_line(const std::string & temp, const std::vector<std::string> & SL) {
	using namespace std;
	for(int i(NUMBER_OF_PREFIXES) ; i < SL.size() ; ++i)
		if(SL[i] == temp)
			return (i - NUMBER_OF_PREFIXES);
	throw(domain_error("ERROR01 in sample_info::find_in_line"));
}


sample_info::sample_info(const std::string & header_line, const std::string & control_filename, const std::string & trio_filename) {
	using namespace std;
	control_index = vector<int>(0);
	trios = vector<vector<int> >(0);
	const vector<string> SL( split(header_line) );
	
	sample_names = vector<string>(SL.begin() + NUMBER_OF_PREFIXES, SL.end());
	
	ifstream fincontrol( control_filename.c_str() );
	if(not fincontrol)
		throw(domain_error("ERROR01 in sample_info::sample_info: control file does not exist"));
	string temp;
	fincontrol >> temp;
	while(fincontrol) {
		control_index.push_back( find_in_line(temp, SL));
		fincontrol >> temp;
	}
	sort(control_index.begin(), control_index.end());
	fincontrol.close();
	
	ifstream fintrio(trio_filename.c_str());
	if(not fintrio)
		throw(domain_error("ERROR02 in sample_info::sample_info: trio file does not exist"));
	getline(fintrio, temp);
	while(fintrio) {
		const vector<string> split_temp( split(temp) );
		if(split_temp.size() != 3)
			throw(domain_error("ERROR03 in sample_info::sample_info: trio file inconsistent format"));
		vector<int> row(0);
		for(int i(0);i<3;++i)
			row.push_back( find_in_line(split_temp[i], SL) );
		trios.push_back( row );
		
		getline(fintrio, temp);
	}
	fintrio.close();
	
	device_control_index = thrust::device_vector<int>(control_index.size());
	thrust::copy(control_index.begin(), control_index.end(), device_control_index.begin());
	
	device_P1_index = thrust::device_vector<int>( trios.size() );
	device_P2_index = thrust::device_vector<int>( trios.size() );
	device_child_index = thrust::device_vector<int>( trios.size() );
	
	for(int i(0);i<trios.size();++i) {
		device_P1_index[i] = trios[i][0];
		device_P2_index[i] = trios[i][1];
		device_child_index[i] = trios[i][2];
	}
}

////////////////

struct EM_worker {
	double rho, lhood, lhood_prev;
	std::vector<double> F;
	size_t iters;
	
	char set_by; // -1 when initialized or reset, 0 when set by GPU computation, 1 when set by CPU computation.
	
	EM_worker();
	
	void reset();
	
	void EM_step    (const AD_info &, const sample_info &);
	void EM_step_CPU(const AD_info &, const sample_info &);
	
	void EM_full    (const AD_info &, const sample_info &);
	void EM_full_CPU(const AD_info &, const sample_info &);
	
};

std::ostream & operator << (std::ostream & out, const EM_worker & O) {
	using namespace std;
	out<< "EM_worker:\n";
	out<< "**********\n";
	if(O.set_by == -1) {
		out<< "initialized or reset\n";
		out<< "--------------------\n";
	}
	else
		if(O.set_by == 0) {
			out<< "set by GPU\n";
			out<< "----------\n";
		}
		else
			if(O.set_by == 1) {
				out<< "set by CPU\n";
				out<< "----------\n";
			}
			else
				throw(domain_error("ERROR in displaying EM_worker"));
	
	out<< "rho = " << O.rho <<"\n";
	out<< "lhood_prev = " << O.lhood_prev <<"\n";
	out<< "lhood = " << O.lhood <<"\n";
	out<< "F: "<< O.F <<"\tsum = "<< SUM(O.F) <<"\n";
	out<< "iters = " << O.iters <<"\n";
	return out;
}

EM_worker::EM_worker() {
	reset();
}

void EM_worker::reset() {
	using namespace std;
	rho = 0.7;
	lhood = 0.0;
	lhood_prev = 0.0;
	F = vector<double>(3, 1.0/3.0);
	iters = 0;
	set_by = -1;
}



struct unary_op {
	double rho;
	double F0, F1, F2;
	double log_rho, log_1_rho;
	double logF0, logF1, logF2;
	unary_op(const double & _rho, const double & _F0, const double & _F1, const double & _F2) {
		rho = _rho;
		F0 = _F0;
		F1 = _F1;
		F2 = _F2;
		log_rho = log(rho);
		log_1_rho = log(1.0 - rho);
		logF0 = log(F0);
		logF1 = log(F1);
		logF2 = log(F2);
	}
	__host__ __device__
	thrust::tuple<double, double, double, double, double, double> operator () (const thrust::tuple<int, int> & DAD) const {
		if( thrust::get<0>(DAD) != 0 or thrust::get<1>(DAD) != 0) {
			double log_w0( thrust::get<0>(DAD)*log_rho + thrust::get<1>(DAD)*log_1_rho + logF0 );
			double log_w1( -( thrust::get<0>(DAD) + thrust::get<1>(DAD) )*LOG_2 + logF1 );
			double log_w2( thrust::get<0>(DAD)*log_1_rho + thrust::get<1>(DAD)*log_rho + logF2 );
			double max_val;
			
			if(log_w0 > log_w1)
				if(log_w0 > log_w2)
					max_val = log_w0;
				else
					max_val = log_w2;
			else
				if(log_w1 > log_w2)
					max_val = log_w1;
				else
					max_val = log_w2;
			
			log_w0 -= max_val;
			log_w1 -= max_val;
			log_w2 -= max_val;
			
			double w0( exp(log_w0) );
			double w1( exp(log_w1) );
			double w2( exp(log_w2) );
			const double den( w0 + w1 + w2 );
			
			const double log_lhood_term( log(den) + max_val );
			
			w0 = w0 / den;
			w1 = w1 / den;
			w2 = w2 / den;
			
			const double T1_term( w0*thrust::get<0>(DAD) + w2*thrust::get<1>(DAD) );
			const double T2_term( w2*thrust::get<0>(DAD) + w0*thrust::get<1>(DAD) );
			
			const thrust::tuple<double, double, double, double, double, double> ans( thrust::make_tuple(T1_term, T2_term, w0, w1, w2, log_lhood_term) );
			return ans;
		}
		else {
			const double w0( 0 );
			const double w1( 0 );
			const double w2( 0 );
			
			const double T1_term( 0 );
			const double T2_term( 0 );
			
			const double log_lhood_term( 0 );
			
			const thrust::tuple<double, double, double, double, double, double> ans( thrust::make_tuple(T1_term, T2_term, w0, w1, w2, log_lhood_term) );
			return ans;
		}
	}
};

struct binary_op {
	__host__ __device__
	thrust::tuple<double, double, double, double, double, double> operator () (const thrust::tuple<double, double, double, double, double, double> & X, const thrust::tuple<double, double, double, double, double, double> & Y) const {
		const double val0( thrust::get<0>(X) + thrust::get<0>(Y) );
		const double val1( thrust::get<1>(X) + thrust::get<1>(Y) );
		const double val2( thrust::get<2>(X) + thrust::get<2>(Y) );
		const double val3( thrust::get<3>(X) + thrust::get<3>(Y) );
		const double val4( thrust::get<4>(X) + thrust::get<4>(Y) );
		const double val5( thrust::get<5>(X) + thrust::get<5>(Y) );
		
		thrust::tuple<double, double, double, double, double, double> ans( thrust::make_tuple(val0, val1, val2, val3, val4, val5) );
		return ans;
	}
};


void EM_worker::EM_step(const AD_info & ADIO, const sample_info & SIO) {
	typedef thrust::device_vector<int>::const_iterator map_iter;
	typedef thrust::device_vector<int>::const_iterator source_iter;
	typedef thrust::permutation_iterator<source_iter, map_iter> perm_iter;
	
	perm_iter beg0   ( ADIO.DAD0.begin(), SIO.device_control_index.begin() );
	perm_iter beg1   ( ADIO.DAD1.begin(), SIO.device_control_index.begin() );
	perm_iter ender0 ( ADIO.DAD0.end()  , SIO.device_control_index.end()   );
	perm_iter ender1 ( ADIO.DAD1.end()  , SIO.device_control_index.end()   );
	
	typedef thrust::tuple<perm_iter, perm_iter> iter_tuple;
	typedef thrust::zip_iterator<iter_tuple> zip_iter;
	
	const zip_iter zip_beg  ( thrust::make_zip_iterator( thrust::make_tuple(beg0  , beg1  ) ) );
	const zip_iter zip_ender( thrust::make_zip_iterator( thrust::make_tuple(ender0, ender1) ) );
	
	const unary_op  UO(rho, F[0], F[1], F[2]);
	const binary_op BO;
	const double lhood_prior_term( (A_prior - 1.0)*log(rho) + (B_prior - 1.0)*log(1.0 - rho) + (D_prior - 1.0)*( log(F[0]) + log(F[1]) + log(F[2]) ) );
	const thrust::tuple<double, double, double, double, double, double> init_val (A_prior - 1.0, B_prior - 1.0, D_prior - 1.0, D_prior - 1.0, D_prior - 1.0, lhood_prior_term);
	const thrust::tuple<double, double, double, double, double, double> TS_VALS( thrust::transform_reduce(zip_beg, zip_ender, UO, init_val, BO) );
	
	const double T1( thrust::get<0>( TS_VALS ) );
	const double T2( thrust::get<1>( TS_VALS ) );
	
	double F0( thrust::get<2>( TS_VALS ) );
	double F1( thrust::get<3>( TS_VALS ) );
	double F2( thrust::get<4>( TS_VALS ) );
	const double den( F0 + F1 + F2 );
	
	F0 = F0 / den;
	F1 = F1 / den;
	F2 = F2 / den;
	
	rho = 1.0/( T2/T1 + 1.0);
	
	F[0] = F0;
	F[1] = F1;
	F[2] = F2;
	
	lhood_prev = lhood;
	lhood = thrust::get<5>(TS_VALS);
}


void EM_worker::EM_step_CPU(const AD_info & ADIO, const sample_info & SIO) {
	using namespace std;
	const thrust::host_vector<int> &AD0( ADIO.AD0 ), &AD1( ADIO.AD1 ), &index( SIO.control_index );
	
	const double lhood_prior_term( (A_prior - 1.0)*log(rho) + (B_prior - 1.0)*log(1.0 - rho) + (D_prior - 1.0)*( log(F[0]) + log(F[1]) + log(F[2]) ) );
	double T1(A_prior - 1.0), T2(B_prior - 1.0), F0(D_prior - 1.0), F1(D_prior - 1.0), F2(D_prior - 1.0);
	lhood_prev = lhood;
	lhood = lhood_prior_term;
	for(int i(0);i<index.size();++i) {
		const int K( index[i] );
		if(AD0[K] != 0 or AD1[K] != 0) {
			double log_w0( AD0[K]*log(rho) + AD1[K]*log(1.0 - rho) + log( F[0] ));
			double log_w1( -(AD0[K] + AD1[K])*LOG_2 + log( F[1] ) );
			double log_w2( AD1[K]*log(rho) + AD0[K]*log(1.0 - rho) + log( F[2] ));
			double max_val;
			
			if(log_w0 > log_w1)
				if(log_w0 > log_w2)
					max_val = log_w0;
				else
					max_val = log_w2;
			else
				if(log_w1 > log_w2)
					max_val = log_w1;
				else
					max_val = log_w2;
			
			log_w0 -= max_val;
			log_w1 -= max_val;
			log_w2 -= max_val;
			
			double w0( exp(log_w0) );
			double w1( exp(log_w1) );
			double w2( exp(log_w2) );
			const double den( w0 + w1 + w2 );
			
			lhood += ( log(den) + max_val );
			
			w0 = w0 / den;
			w1 = w1 / den;
			w2 = w2 / den;
			
			const double T1_term( w0*AD0[K] + w2*AD1[K] );
			const double T2_term( w2*AD0[K] + w0*AD1[K] );
			
			T1 += T1_term;
			T2 += T2_term;
			F0 += w0;
			F1 += w1;
			F2 += w2;
		}
		else {
			T1 += 0.0;
			T2 += 0.0;
			F0 += 0.0;
			F1 += 0.0;
			F2 += 0.0;
			lhood += 0.0;
		}
	}
	const double den( F0 + F1 + F2 );
	
	F0 = F0 / den;
	F1 = F1 / den;
	F2 = F2 / den;
	
	rho = 1.0/( T2/T1 + 1.0);
	F[0] = F0;
	F[1] = F1;
	F[2] = F2;
}


void EM_worker::EM_full(const AD_info & ADIO, const sample_info & SIO) {
	reset();
	EM_step(ADIO, SIO);
	EM_step(ADIO, SIO);
	iters = 2;
	while(abs(lhood - lhood_prev) > THRESHOLD) {
		EM_step(ADIO, SIO);
		iters += 1;
		if(iters > MAX_ITER)
			break;
	}
	set_by = 0;
}

void EM_worker::EM_full_CPU(const AD_info & ADIO, const sample_info & SIO) {
	reset();
	EM_step_CPU(ADIO, SIO);
	EM_step_CPU(ADIO, SIO);
	iters = 2;
	while(abs(lhood - lhood_prev) > THRESHOLD) {
		EM_step_CPU(ADIO, SIO);
		iters += 1;
		if(iters > MAX_ITER)
			break;
	}
	set_by = 1;
}


//////////////////

struct conditional_probty_table {
	std::vector<double> CPT; // conditional_probty_table
	thrust::device_vector<double> device_CPT; // GPU memory conditional_probty_table
	
	std::vector<double> LOG_CPT; // log of conditional_probty_table
	thrust::device_vector<double> device_LOG_CPT; // GPU memory log of conditional_probty_table
	
	conditional_probty_table();
};

conditional_probty_table::conditional_probty_table() {
		double temp[] = {1,0,0,      0.5,0.5,0,      0,1,0,      0.5,0.5,0,      0.25,0.5,0.25,      0,0.5,0.5,      0,1,0,      0,0.5,0.5,      0,0,1};
		
		CPT = std::vector<double>(temp, temp + 27);
		CPT += MUTATION_RATE;
		LOG_CPT = LOG_1D(CPT);
		
		device_CPT = thrust::device_vector<double>( CPT );
		device_LOG_CPT = thrust::device_vector<double>( LOG_CPT );
}

std::ostream & operator << (std::ostream & out, const conditional_probty_table & O) {
	out << "conditional_probty_table:\n";
	out << "*************************\n";
	out << "CPT:\n";
	out << O.CPT <<"\t"<< ":length = "<< O.CPT.size() <<"\n";
	out << "LOG_CPT:\n";
	out << O.LOG_CPT <<"\t"<< ":length = "<< O.LOG_CPT.size() <<"\n";
	out << "device_CPT:\n";
	out << O.device_CPT <<"\t"<< ":length = "<< O.device_CPT.size() <<"\n";
	out << "device_LOG_CPT:\n";
	out << O.device_LOG_CPT <<"\t"<< ":length = "<< O.device_LOG_CPT.size() <<"\n";
	return out;
}

//////////////////

struct caller_unary_op {
	double CPTL000,  CPTL001,  CPTL002,  CPTL010,  CPTL011,  CPTL012,  CPTL020,  CPTL021,  CPTL022,  CPTL100,  CPTL101,  CPTL102,  CPTL110,  CPTL111,  CPTL112,  CPTL120,  CPTL121,  CPTL122,  CPTL200,  CPTL201,  CPTL202,  CPTL210,  CPTL211,  CPTL212,  CPTL220,  CPTL221,  CPTL222;
	double rho;
	double F0, F1, F2;
	
	caller_unary_op(const conditional_probty_table & O, const double & _rho, const std::vector<double> & _F);
	__host__ __device__
	double operator () (const thrust::tuple<int, int, int, int, int, int> & X) const;
};

caller_unary_op::caller_unary_op(const conditional_probty_table & O, const double & _rho, const std::vector<double> & _F) {
	if(O.LOG_CPT.size() != 27)
		throw(std::domain_error("ERROR in caller_unary_op::caller_unary_op"));
	
	CPTL000 = O.LOG_CPT[0];
	CPTL001 = O.LOG_CPT[1];
	CPTL002 = O.LOG_CPT[2];
	CPTL010 = O.LOG_CPT[3];
	CPTL011 = O.LOG_CPT[4];
	CPTL012 = O.LOG_CPT[5];
	CPTL020 = O.LOG_CPT[6];
	CPTL021 = O.LOG_CPT[7];
	CPTL022 = O.LOG_CPT[8];
	CPTL100 = O.LOG_CPT[9];
	CPTL101 = O.LOG_CPT[10];
	CPTL102 = O.LOG_CPT[11];
	CPTL110 = O.LOG_CPT[12];
	CPTL111 = O.LOG_CPT[13];
	CPTL112 = O.LOG_CPT[14];
	CPTL120 = O.LOG_CPT[15];
	CPTL121 = O.LOG_CPT[16];
	CPTL122 = O.LOG_CPT[17];
	CPTL200 = O.LOG_CPT[18];
	CPTL201 = O.LOG_CPT[19];
	CPTL202 = O.LOG_CPT[20];
	CPTL210 = O.LOG_CPT[21];
	CPTL211 = O.LOG_CPT[22];
	CPTL212 = O.LOG_CPT[23];
	CPTL220 = O.LOG_CPT[24];
	CPTL221 = O.LOG_CPT[25];
	CPTL222 = O.LOG_CPT[26];
	
	rho = _rho;
	F0 = _F[0];
	F1 = _F[1];
	F2 = _F[2];
}

double caller_unary_op::operator () (const thrust::tuple<int, int, int, int, int, int> & X) const {
	// parent1 ADs, parent2 ADs, child ADs in that order
	const int P1_AD0( thrust::get<0>(X) );
	const int P1_AD1( thrust::get<1>(X) );
	
	const int P2_AD0( thrust::get<2>(X) );
	const int P2_AD1( thrust::get<3>(X) );
	
	const int child_AD0( thrust::get<4>(X) );
	const int child_AD1( thrust::get<5>(X) );
	
	thrust::tuple<double, double, double> P1L, P2L, childL;
	{
		thrust::get<0>(P1L) = P1_AD0*log(rho) + P1_AD1*log(1.0 - rho);
		thrust::get<1>(P1L) = -(P1_AD0 + P1_AD1)*LOG_2;
		thrust::get<2>(P1L) = P1_AD1*log(rho) + P1_AD0*log(1.0 - rho);
		
		thrust::get<0>(P2L) = P2_AD0*log(rho) + P2_AD1*log(1.0 - rho);
		thrust::get<1>(P2L) = -(P2_AD0 + P2_AD1)*LOG_2;
		thrust::get<2>(P2L) = P2_AD1*log(rho) + P2_AD0*log(1.0 - rho);
		
		thrust::get<0>(childL) = child_AD0*log(rho) + child_AD1*log(1.0 - rho);
		thrust::get<1>(childL) = -(child_AD0 + child_AD1)*LOG_2;
		thrust::get<2>(childL) = child_AD1*log(rho) + child_AD0*log(1.0 - rho);
	}
	
	double PP000,  PP001,  PP002,  PP010,  PP011,  PP012,  PP020,  PP021,  PP022,  PP100,  PP101,  PP102,  PP110,  PP111,  PP112,  PP120,  PP121,  PP122,  PP200,  PP201,  PP202,  PP210,  PP211,  PP212,  PP220,  PP221,  PP222;
	
	PP000 = CPTL000 + thrust::get<0>(P1L) + thrust::get<0>(P2L) + thrust::get<0>(childL);
	PP001 = CPTL001 + thrust::get<0>(P1L) + thrust::get<0>(P2L) + thrust::get<1>(childL);
	PP002 = CPTL002 + thrust::get<0>(P1L) + thrust::get<0>(P2L) + thrust::get<2>(childL);
	PP010 = CPTL010 + thrust::get<0>(P1L) + thrust::get<1>(P2L) + thrust::get<0>(childL);
	PP011 = CPTL011 + thrust::get<0>(P1L) + thrust::get<1>(P2L) + thrust::get<1>(childL);
	PP012 = CPTL012 + thrust::get<0>(P1L) + thrust::get<1>(P2L) + thrust::get<2>(childL);
	PP020 = CPTL020 + thrust::get<0>(P1L) + thrust::get<2>(P2L) + thrust::get<0>(childL);
	PP021 = CPTL021 + thrust::get<0>(P1L) + thrust::get<2>(P2L) + thrust::get<1>(childL);
	PP022 = CPTL022 + thrust::get<0>(P1L) + thrust::get<2>(P2L) + thrust::get<2>(childL);
	PP100 = CPTL100 + thrust::get<1>(P1L) + thrust::get<0>(P2L) + thrust::get<0>(childL);
	PP101 = CPTL101 + thrust::get<1>(P1L) + thrust::get<0>(P2L) + thrust::get<1>(childL);
	PP102 = CPTL102 + thrust::get<1>(P1L) + thrust::get<0>(P2L) + thrust::get<2>(childL);
	PP110 = CPTL110 + thrust::get<1>(P1L) + thrust::get<1>(P2L) + thrust::get<0>(childL);
	PP111 = CPTL111 + thrust::get<1>(P1L) + thrust::get<1>(P2L) + thrust::get<1>(childL);
	PP112 = CPTL112 + thrust::get<1>(P1L) + thrust::get<1>(P2L) + thrust::get<2>(childL);
	PP120 = CPTL120 + thrust::get<1>(P1L) + thrust::get<2>(P2L) + thrust::get<0>(childL);
	PP121 = CPTL121 + thrust::get<1>(P1L) + thrust::get<2>(P2L) + thrust::get<1>(childL);
	PP122 = CPTL122 + thrust::get<1>(P1L) + thrust::get<2>(P2L) + thrust::get<2>(childL);
	PP200 = CPTL200 + thrust::get<2>(P1L) + thrust::get<0>(P2L) + thrust::get<0>(childL);
	PP201 = CPTL201 + thrust::get<2>(P1L) + thrust::get<0>(P2L) + thrust::get<1>(childL);
	PP202 = CPTL202 + thrust::get<2>(P1L) + thrust::get<0>(P2L) + thrust::get<2>(childL);
	PP210 = CPTL210 + thrust::get<2>(P1L) + thrust::get<1>(P2L) + thrust::get<0>(childL);
	PP211 = CPTL211 + thrust::get<2>(P1L) + thrust::get<1>(P2L) + thrust::get<1>(childL);
	PP212 = CPTL212 + thrust::get<2>(P1L) + thrust::get<1>(P2L) + thrust::get<2>(childL);
	PP220 = CPTL220 + thrust::get<2>(P1L) + thrust::get<2>(P2L) + thrust::get<0>(childL);
	PP221 = CPTL221 + thrust::get<2>(P1L) + thrust::get<2>(P2L) + thrust::get<1>(childL);
	PP222 = CPTL222 + thrust::get<2>(P1L) + thrust::get<2>(P2L) + thrust::get<2>(childL);
	
	double maxval( PP000 );
	if(maxval < PP000)
		maxval = PP000;
	if(maxval < PP001)
		maxval = PP001;
	if(maxval < PP002)
		maxval = PP002;
	if(maxval < PP010)
		maxval = PP010;
	if(maxval < PP011)
		maxval = PP011;
	if(maxval < PP012)
		maxval = PP012;
	if(maxval < PP020)
		maxval = PP020;
	if(maxval < PP021)
		maxval = PP021;
	if(maxval < PP022)
		maxval = PP022;
	if(maxval < PP100)
		maxval = PP100;
	if(maxval < PP101)
		maxval = PP101;
	if(maxval < PP102)
		maxval = PP102;
	if(maxval < PP110)
		maxval = PP110;
	if(maxval < PP111)
		maxval = PP111;
	if(maxval < PP112)
		maxval = PP112;
	if(maxval < PP120)
		maxval = PP120;
	if(maxval < PP121)
		maxval = PP121;
	if(maxval < PP122)
		maxval = PP122;
	if(maxval < PP200)
		maxval = PP200;
	if(maxval < PP201)
		maxval = PP201;
	if(maxval < PP202)
		maxval = PP202;
	if(maxval < PP210)
		maxval = PP210;
	if(maxval < PP211)
		maxval = PP211;
	if(maxval < PP212)
		maxval = PP212;
	if(maxval < PP220)
		maxval = PP220;
	if(maxval < PP221)
		maxval = PP221;
	if(maxval < PP222)
		maxval = PP222;
	
	PP000 -= maxval;
	PP001 -= maxval;
	PP002 -= maxval;
	PP010 -= maxval;
	PP011 -= maxval;
	PP012 -= maxval;
	PP020 -= maxval;
	PP021 -= maxval;
	PP022 -= maxval;
	PP100 -= maxval;
	PP101 -= maxval;
	PP102 -= maxval;
	PP110 -= maxval;
	PP111 -= maxval;
	PP112 -= maxval;
	PP120 -= maxval;
	PP121 -= maxval;
	PP122 -= maxval;
	PP200 -= maxval;
	PP201 -= maxval;
	PP202 -= maxval;
	PP210 -= maxval;
	PP211 -= maxval;
	PP212 -= maxval;
	PP220 -= maxval;
	PP221 -= maxval;
	PP222 -= maxval;
	
	double num( exp(PP001) );
	double den(0);
	
	den += exp(PP000);
	den += exp(PP001);
	den += exp(PP002);
	den += exp(PP010);
	den += exp(PP011);
	den += exp(PP012);
	den += exp(PP020);
	den += exp(PP021);
	den += exp(PP022);
	den += exp(PP100);
	den += exp(PP101);
	den += exp(PP102);
	den += exp(PP110);
	den += exp(PP111);
	den += exp(PP112);
	den += exp(PP120);
	den += exp(PP121);
	den += exp(PP122);
	den += exp(PP200);
	den += exp(PP201);
	den += exp(PP202);
	den += exp(PP210);
	den += exp(PP211);
	den += exp(PP212);
	den += exp(PP220);
	den += exp(PP221);
	den += exp(PP222);
	
	return (num/den);
}

//////////////////

struct caller {
	char set_by; // -1 when initialized or reset, 0 when set by GPU computation, 1 when set by CPU computation.
	thrust::host_vector<double> call_probty;
	thrust::device_vector<double> device_call_probty;
	
	caller(const sample_info &);
	
	void reset();
	
	void worker    (const AD_info & ADIO, const EM_worker & EMWO, const sample_info & SIO, const conditional_probty_table & CPTO);
	void worker_CPU(const AD_info & ADIO, const EM_worker & EMWO, const sample_info & SIO, const conditional_probty_table & CPTO);
	
	private:
	double worker_CPU_aux(const std::vector<std::vector<int> > &, const EM_worker &, const conditional_probty_table &);
};

void caller::reset() {
	for(size_t i(0);i<call_probty.size();++i)
		call_probty[i] = 0;
	
	thrust::copy(call_probty.begin(), call_probty.end(), device_call_probty.begin());
	set_by = -1;
}

std::ostream & operator << (std::ostream & out, const caller & O) {
	using namespace std;
	out << "caller:\n";
	out << "*******\n";
	if(O.set_by == -1) {
		out<< "initialized or reset\n";
		out<< "--------------------\n";
	}
	else
		if(O.set_by == 0) {
			out<< "set by GPU\n";
			out<< "----------\n";
		}
		else
			if(O.set_by == 1) {
				out<< "set by CPU\n";
				out<< "----------\n";
			}
			else
				throw(domain_error("ERROR in displaying caller"));
	
	out << "call_probty:\n";
	out << O.call_probty <<"\t"<< ":length = " << O.call_probty.size() <<"\n";
	out << "device_call_probty:\n";
	out << O.device_call_probty <<"\t"<< ":length = "<< O.device_call_probty.size() <<"\n";
	return out;
}

caller::caller(const sample_info & SIO) {
	call_probty = thrust::host_vector<double>( SIO.trios.size() );
	device_call_probty = thrust::device_vector<double>( SIO.trios.size() );
	set_by = -1;
}


void caller::worker(const AD_info & ADIO, const EM_worker & EMWO, const sample_info & SIO, const conditional_probty_table & CPTO) {
	typedef thrust::device_vector<int>::const_iterator map_iter;
	typedef thrust::device_vector<int>::const_iterator source_iter;
	typedef thrust::permutation_iterator<source_iter, map_iter> perm_iter;
	
	perm_iter P1_beg0     ( ADIO.DAD0.begin() , SIO.device_P1_index.begin()    );
	perm_iter P1_beg1     ( ADIO.DAD1.begin() , SIO.device_P1_index.begin()    );
	
	perm_iter P1_ender0   ( ADIO.DAD0.end()   , SIO.device_P1_index.end()      );
	perm_iter P1_ender1   ( ADIO.DAD1.end()   , SIO.device_P1_index.end()      );
	
	perm_iter P2_beg0     ( ADIO.DAD0.begin() , SIO.device_P2_index.begin()    );
	perm_iter P2_beg1     ( ADIO.DAD1.begin() , SIO.device_P2_index.begin()    );
	
	perm_iter P2_ender0   ( ADIO.DAD0.end()   , SIO.device_P2_index.end()      );
	perm_iter P2_ender1   ( ADIO.DAD1.end()   , SIO.device_P2_index.end()      );
	
	perm_iter child_beg0   ( ADIO.DAD0.begin(), SIO.device_child_index.begin() );
	perm_iter child_beg1   ( ADIO.DAD1.begin(), SIO.device_child_index.begin() );
	
	perm_iter child_ender0 ( ADIO.DAD0.end()  , SIO.device_child_index.end()   );
	perm_iter child_ender1 ( ADIO.DAD1.end()  , SIO.device_child_index.end()   );
	
	typedef thrust::tuple<perm_iter, perm_iter, perm_iter, perm_iter, perm_iter, perm_iter> iter_tuple;
	typedef thrust::zip_iterator<iter_tuple> zip_iter;
	
	zip_iter zip_beg  ( thrust::make_zip_iterator( thrust::make_tuple(P1_beg0  , P1_beg1  , P2_beg0  , P2_beg1  , child_beg0  , child_beg1  ) ) );
	zip_iter zip_ender( thrust::make_zip_iterator( thrust::make_tuple(P1_ender0, P1_ender1, P2_ender0, P2_ender1, child_ender0, child_ender1) ) );
	
	caller_unary_op CUOO(CPTO, EMWO.rho, EMWO.F);
	
	thrust::transform(zip_beg, zip_ender, device_call_probty.begin(), CUOO);
	
	thrust::copy(device_call_probty.begin(), device_call_probty.end(), call_probty.begin());
	
	set_by = 0;
}


double caller::worker_CPU_aux(const std::vector<std::vector<int> > & trio_ADs, const EM_worker & EMWO, const conditional_probty_table & CPTO) {
	using namespace std;
	const double rho( EMWO.rho);
	const vector<double> F( EMWO.F );
	vector<vector<double> > LLs(3, vector<double>(3));
	
	for(int i(0);i<3;++i) {
		vector<double> & row( LLs[i] );
		const vector<int> & AD( trio_ADs[i] );
		if(AD.size() != 2)
			throw(domain_error("ERROR in worker_CPU_aux"));
		row[0] = AD[0]*log(rho) + AD[1]*log(1.0 - rho);
		row[1] = -(AD[0] + AD[1])*LOG_2;
		row[2] = AD[1]*log(rho) + AD[0]*log(1.0 - rho);
	}
	
	
	vector<double> combined(27);
	for(int i(0);i<3;++i) {
		for(int j(0);j<3;++j) {
			for(int k(0);k<3;++k) {
				const int I(9*i + 3*j + k);
				combined[I] = CPTO.LOG_CPT[I] + LLs[0][i] + LLs[1][j] + LLs[2][k];
			}
		}
	}
	
	const double maxval( MAX(combined) );
	combined += (-maxval);
	EXP_1D_INPLACE(combined);
	const double den( SUM(combined) );
	
	return ( combined[1]/den );
}

void caller::worker_CPU(const AD_info & ADIO, const EM_worker & EMWO, const sample_info & SIO, const conditional_probty_table & CPTO) {
	using namespace std;
	if(SIO.trios.size() != call_probty.size())
		throw(domain_error("ERROR in caller::worker_CPU"));
	
	for(int i(0);i<SIO.trios.size();++i) {
		const vector<int> Index( SIO.trios[i] );
		
		vector<int> P1_AD(2);
		P1_AD[0] = ADIO.AD0[ Index[0] ];
		P1_AD[1] = ADIO.AD1[ Index[0] ];
		vector<int> P2_AD(2);
		P2_AD[0] = ADIO.AD0[ Index[1] ];
		P2_AD[1] = ADIO.AD1[ Index[1] ];
		vector<int> child_AD(2);
		child_AD[0] = ADIO.AD0[ Index[2] ];
		child_AD[1] = ADIO.AD1[ Index[2] ];
		
		vector<vector<int> > trio_ADs(3);
		trio_ADs[0] = P1_AD;
		trio_ADs[1] = P2_AD;
		trio_ADs[2] = child_AD;	
		
		const double ans( worker_CPU_aux(trio_ADs, EMWO, CPTO) );
		call_probty[i] = ans;
	}
	
	thrust::copy(call_probty.begin(), call_probty.end(), device_call_probty.begin());
	
	set_by = 1;
}



//////////////////


#endif

















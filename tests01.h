#ifndef GAURD_tests01
#define GAURD_tests01



void tester_01() {
	using namespace std;
	cout.precision(15);
	ifstream fin("/home/anwoy/useful_data/vcf_data/bgm0042.vep.hg19.vcf");
	if(not fin)
		throw(domain_error("ERROR: file not present"));
	
	string line;
	getline(fin, line);
	while(true) {
		if(line.size() >= 6) {
			const string temp(line.begin(), line.begin() + 6);
			if(temp == "#CHROM")
				break;
		}
		getline(fin, line);
	}
	
	const sample_info SIO(line, "controls.txt", "trios01.txt");
	cout<< SIO <<"\n";
	const initialization_info IIO(line);
	cout<< IIO <<"\n";
	const conditional_probty_table CPTO;
	cout<< CPTO <<"\n";
	
	AD_info ADIO(IIO);
	size_t counter(0);
	while(ADIO.read_once(fin, IIO)) {
		counter += 1;
		if(counter == 100)
			break;
	}
	cout<< ADIO <<"\n";
	
	EM_worker EMWO;
	cout<< EMWO <<"\n";
	const clock_t T1( clock() );
	EMWO.EM_full(ADIO, SIO); // GPU calculation
	const clock_t T2( clock() );
	cout<< "GPU time = "<< double(T2 - T1)/CLOCKS_PER_SEC <<" seconds\n";
	cout<< EMWO <<"\n";
	EMWO.reset();
	cout<< EMWO <<"\n";
	const clock_t T3( clock() );
	EMWO.EM_full_CPU(ADIO, SIO); // CPU calculation
	const clock_t T4( clock() );
	cout<< "CPU time = "<< double(T4 - T3)/CLOCKS_PER_SEC <<" seconds\n";
	cout<< EMWO <<"\n";
	
	caller CO(SIO);
	cout<< CO <<"\n";
	const clock_t T5( clock() );
	CO.worker(ADIO, EMWO, SIO, CPTO);
	const clock_t T6( clock() );
	cout<< CO <<"\n";
	CO.reset();
	cout<< CO <<"\n";
	const clock_t T7( clock() );
	CO.worker_CPU(ADIO, EMWO, SIO, CPTO);
	const clock_t T8( clock() );
	cout<< CO <<"\n";
	
	cout<< "GPU time for caller = " << double(T6 - T5)/CLOCKS_PER_SEC <<"\n";
	cout<< "CPU time for caller = " << double(T8 - T7)/CLOCKS_PER_SEC <<"\n";
}

void tester_02() {
	using namespace std;
	cout.precision(15);
	ifstream fin("/home/anwoy/useful_data/vcf_data/bgm0042.vep.hg19.vcf");
	if(not fin)
		throw(domain_error("ERROR: file not present"));
	
	string line;
	getline(fin, line);
	while(true) {
		if(line.size() >= 6) {
			const string temp(line.begin(), line.begin() + 6);
			if(temp == "#CHROM")
				break;
		}
		getline(fin, line);
	}
	
	const sample_info SIO(line, "controls.txt", "trios.txt");
	const initialization_info IIO(line);
	const conditional_probty_table CPTO;
	EM_worker EMWO;
	caller CO(SIO);
	AD_info ADIO(IIO);
	
	int counter(0);
	const clock_t T1( clock() );
	while(ADIO.read_once(fin, IIO)) {
		counter += 1;
		EMWO.EM_full_CPU(ADIO, SIO); // CPU calculation
		CO.worker_CPU(ADIO, EMWO, SIO, CPTO); // CPU calculation
		
		if(counter % 1000000 == 0) {
			const clock_t T2( clock() );
			cout<< counter <<" "<< double(T2 - T1)/CLOCKS_PER_SEC <<"\n"; ;
			cout.flush();
		}
		
		bool mark(0);
		for(int i(0);i<CO.call_probty.size();++i)
			if(CO.call_probty[i] > 0.9)
				mark = 1;
		if(mark)
			break;
	}
	cout<< SIO <<"\n";
	cout<< IIO <<"\n";
	cout<< CPTO <<"\n";
	cout<< ADIO <<"\n";
	
	cout<< EMWO <<"\n";
	cout<< CO <<"\n";
	EMWO.reset();
	CO.reset();
	cout<< EMWO <<"\n";
	cout<< CO <<"\n";
	EMWO.EM_full(ADIO, SIO); // GPU calculation
	CO.worker(ADIO, EMWO, SIO, CPTO); // GPU calculation
	cout<< EMWO <<"\n";
	cout<< CO <<"\n";
}

void master_CPU(const std::string & vcf_file, const std::string & controls_file, const std::string & trios_file, const std::string & outfile,const double THRESH) {
	// computes the denovo calls using CPU
	using namespace std;
	cout<<"computing using CPU:\n";
	
	cout<<"vcf_file = "<< vcf_file <<"\n";
	cout<<"controls_file = "<< controls_file <<"\n";
	cout<<"trios_file = "<< trios_file <<"\n";
	cout<<"probability threshold = "<< THRESH <<"\n";
	
	cout.precision(15);
	ifstream fin(vcf_file.c_str());
	if(not fin)
		throw(domain_error("ERROR: file not present"));
	
	string line;
	getline(fin, line);
	while(true) {
		if(line.size() >= 6) {
			const string temp(line.begin(), line.begin() + 6);
			if(temp == "#CHROM")
				break;
		}
		getline(fin, line);
	}
	
	const sample_info SIO(line, controls_file, trios_file);
	const initialization_info IIO(line);
	const conditional_probty_table CPTO;
	EM_worker EMWO;
	caller CO(SIO);
	AD_info ADIO(IIO);
	
	int counter(0);
	const clock_t T1( clock() );
	ofstream fout(outfile.c_str());
	fout.precision(10);
	while(ADIO.read_once(fin, IIO)) {
		counter += 1;
		if(ADIO.CHROM != 24 and ADIO.CHROM != 0) {
			EMWO.EM_full_CPU(ADIO, SIO); // CPU calculation
			CO.worker_CPU(ADIO, EMWO, SIO, CPTO); // CPU calculation
			
			bool mark(0);
			for(int i(0);i<CO.call_probty.size();++i) {
				if(CO.call_probty[i] > THRESH) {
					if(mark == 0) {
						fout<< ADIO.integer_to_chromosome() <<"\t"<< ADIO.POS <<"\t";
						mark = 1;
					}
					const vector<int> & row( SIO.trios[i] );
					fout<< row[0] <<":"<< SIO.sample_names[ row[0] ] <<",";
					fout<< row[1] <<":"<< SIO.sample_names[ row[1] ] <<",";
					fout<< row[2] <<":"<< SIO.sample_names[ row[2] ] <<",";
					fout<< "PP=" << CO.call_probty[i] <<"\t";
				}
			}
			if(mark == 1)
				fout<<"\n";
		}
		
		if(counter % 1000000 == 0) {
			const clock_t T2( clock() );
			cout<< counter <<" lines\t"<< double(T2 - T1)/CLOCKS_PER_SEC <<" seconds\n";
		}
	}
	
	fin.close();
	fout.close();
}

void master(const std::string & vcf_file, const std::string & controls_file, const std::string & trios_file, const std::string & outfile,const double THRESH) {
	// computes the denovo calls using GPU
	using namespace std;
	cout<<"computing using GPU:\n";
	
	cout<<"vcf_file = "<< vcf_file <<"\n";
	cout<<"controls_file = "<< controls_file <<"\n";
	cout<<"trios_file = "<< trios_file <<"\n";
	cout<<"probability threshold = "<< THRESH <<"\n";
	
	cout.precision(15);
	ifstream fin(vcf_file.c_str());
	if(not fin)
		throw(domain_error("ERROR: file not present"));
	
	string line;
	getline(fin, line);
	while(true) {
		if(line.size() >= 6) {
			const string temp(line.begin(), line.begin() + 6);
			if(temp == "#CHROM")
				break;
		}
		getline(fin, line);
	}
	
	const sample_info SIO(line, controls_file, trios_file);
	const initialization_info IIO(line);
	const conditional_probty_table CPTO;
	EM_worker EMWO;
	caller CO(SIO);
	AD_info ADIO(IIO);
	
	int counter(0);
	const clock_t T1( clock() );
	ofstream fout(outfile.c_str());
	fout.precision(10);
	while(ADIO.read_once(fin, IIO)) {
		counter += 1;
		if(ADIO.CHROM != 24 and ADIO.CHROM != 0) {
			EMWO.EM_full(ADIO, SIO); // GPU calculation
			CO.worker(ADIO, EMWO, SIO, CPTO); // GPU calculation
			
			bool mark(0);
			for(int i(0);i<CO.call_probty.size();++i) {
				if(CO.call_probty[i] > THRESH) {
					if(mark == 0) {
						fout<< ADIO.integer_to_chromosome() <<"\t"<< ADIO.POS <<"\t";
						mark = 1;
					}
					const vector<int> & row( SIO.trios[i] );
					fout<< row[0] <<":"<< SIO.sample_names[ row[0] ] <<",";
					fout<< row[1] <<":"<< SIO.sample_names[ row[1] ] <<",";
					fout<< row[2] <<":"<< SIO.sample_names[ row[2] ] <<",";
					fout<< "PP=" << CO.call_probty[i] <<"\t";
				}
			}
			if(mark == 1)
				fout<<"\n";
		}
		
		if(counter % 100000 == 0) {
			const clock_t T2( clock() );
			cout<< counter <<" lines\t"<< double(T2 - T1)/CLOCKS_PER_SEC <<" seconds\n";
		}
	}
	
	fin.close();
	fout.close();
}


void tester_03() {
	using namespace std;
	const string vcf_file("/home/anwoy/useful_data/vcf_data/bgm0042.vep.hg19.vcf");
	const string controls_file("controls.txt");
	const string trios_file("trios02.txt");
	const string outfile_CPU("outfile_CPU.txt");
	const string outfile_GPU("outfile_GPU.txt");
	const double THRESH(0.97);
	
	//master_CPU(vcf_file, controls_file, trios_file, outfile_CPU, THRESH); // CPU computation
	master(vcf_file, controls_file, trios_file, outfile_GPU, THRESH); // GPU computation
}


#endif


















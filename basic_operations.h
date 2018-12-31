#ifndef GAURD_basic_operations
#define GAURD_basic_operations

//using namespace std;

template<class T>
std::ostream& operator<<(std::ostream &out,const std::vector<T> &V){
	using namespace std;
	if(V.size() > 0){
		for(size_t i(0);i<V.size()-1;++i)
			out<<V[i]<<",";
		out<<V[V.size()-1];
	}
	return out;
}

template<class T>
std::ostream& operator<<(std::ostream &out,const std::vector<std::vector<T> > &V){
	using namespace std;
	for(int i(0);i<V.size();++i)
		out<<V[i]<<"\n";
	return out;
}

template<class T>
std::ostream & operator << (std::ostream & out, const thrust::device_vector<T> & V) {
	for(int i(0);i<V.size()-1;++i)
		out<<V[i]<<",";
	out<<V[V.size()-1];
	return out;
}

template<class T>
std::ostream & operator << (std::ostream & out, const thrust::host_vector<T> & V) {
	for(int i(0);i<V.size()-1;++i)
		out<<V[i]<<",";
	out<<V[V.size()-1];
	return out;
}


std::vector<std::string> split(const std::string& line){
    std::vector<std::string> split_string;
    std::string::const_iterator i(line.begin()),j(line.begin());
    while(i!=line.end() && j!=line.end()){
        while(isspace(*i) && i!=line.end())
            ++i;
        j=i;
        while(!isspace(*j) && j!=line.end())
            ++j;
        std::string word(i,j);
        if(word.length()>0)
            split_string.push_back(word);
        i=j;
    }
    return split_string;
}

std::vector<std::string> split(const std::string & line,const std::string & S) {
	using namespace std;
	vector<string> ans;
	const size_t L( S.size() );
	string::const_iterator i(line.begin()),j1(line.begin()),j2(line.begin() + L);
	while(i != line.end()){
		{
			const string temp(i,line.end());
			if(temp.size() < L){
				ans.push_back(temp);
				return ans;
			}
		}
		
		j1 = i;
		j2 = j1+L;
		string val(j1,j2);
		while(val != S && j2 != line.end()){
			++j1;
			j2 = j1+L;
			val = string(j1,j2);
		}
		
		if(val == S){
			string temp(i,j1);
			ans.push_back(temp);
		}
		else{
			string temp(i,line.end());
			ans.push_back(temp);
			return ans;
		}
		
		if(j2 == line.end()){
			ans.push_back("");
			return ans;
		}
		
		i = j2;
	}
	
	return ans;
}

std::string myitoa(int CHROM) {
	using namespace std;
	string sign("");
	if(CHROM < 0) {
		sign = "-";
		CHROM = -CHROM;
	}
	string val("");
	while(CHROM > 0) {
		const int rem(CHROM % 10);
		const string temp(1, '0'+rem);
		val = temp + val;
		CHROM /= 10;
	}
	val = sign + val;
	return val;
}






#endif































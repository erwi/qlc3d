#include <reader.h>

Reader::Reader()
{
    //file = new std::fstream();
    isIgnoreCase = false;
    allowSpaces = false;
    //ctor
}

Reader::~Reader()
{
    //dtor
}

bool Reader::openFile(const std::string& filename){

    file.open( filename.c_str(), std::fstream::in | std::fstream::out );

    return file.is_open();
}
// end openFile
void Reader::closeFile(){
    if (file.is_open() )
	file.close();

}
// end closeFile

std::string Reader::getErrorString(int e){
    using namespace std;
    switch(e){
        case READER_SUCCESS:{
                return std::string("READER_SUCCESS");
                break;
            }
        case READER_NOT_FOUND:{
            return string("READER_NOT_FOUND");
            break;
        }
        case READER_BAD_VALUE:{
            return string("READER_BAD_VALUE");
            break;
        }
        case READER_BAD_FORMAT:{
            return string("READER_BAD_FORMAT");
            break;
        }
        case READER_BAD_FILE:{
			return string("READER_BAD_FILE");
			break;
		}
        default:{
            return string("UNDEFINED ERROR");
            break;
            }
        }

}
//end getErrorString

int Reader::findStringLine(std::string str,
			   unsigned int offset  ){
	file.seekg(0);
	gotoLine( offset );
	if (file.good() ){
		std::string line;
		stringToLowercase(str);
		//std::cout << "searching for " << str << std::endl;
		int lc = 0; //line counter
		while( !file.eof() ){
			getline(file , line);
			removeBlanksAtBeginning(line);
			removeBlanksAtEnd(line);
			stringToLowercase(line);
			//std::cout << "comparing with " << line << std::endl;
			if (!line.compare(str)){ // string found
			    //std::cout << str << "on line number" << lc << std::endl;
				return lc;
			}// end if string found
			lc++;
		}// end while !eof
	}// file good?
	else{
		return READER_BAD_FILE;
	}// end if file good
	return READER_BAD_FILE;
}

int Reader::findVariableLine(const std::string& var,   // returns line next number where variable var is assigned
                     unsigned int offset){ // if offset ==0 , line number is w.r.t. beginning of file
	using std::cout;
	using std::endl;
	file.seekg(0);
	gotoLine( offset );

	if (!file.good() )
		return READER_BAD_FILE;

	int lc = 0; // line count
	std::string line;
	std::string svar = var;

	if (isIgnoreCase)			// if no case sensitivity, make lower case
		stringToLowercase(svar);

	while( ! file.eof() ){
		getline(file, line );
		removeBlanks(line);
		if (isIgnoreCase)
			stringToLowercase(line);

		if ( line.find(svar) == 0 ){ // if line starts with variable

			// further tests are needed. For a valid assingment, next character must be '='
			size_t len = svar.length();
			if (line.find_first_of('=') == len ){
				return lc;
			}

		}
		lc++;
	}// end of line
	file.seekg(0);
	file.clear();
	return READER_NOT_FOUND;

}
int Reader::readString(std::string var, std::string& val){
    // std::cout << "searching for :" << var << std::endl;
    file.seekg(0);
    if (!file.good() )
        return READER_BAD_FILE;

	int ln = this->findVariableLine(var);


	if (ln >= 0 ){ // if found
		this->gotoLine( (unsigned int) ln );
		std::string line;
                getline(file, line);

                
                val = this->extractRHS(line);
                //printf(" RHS =%s. \n", val.c_str() );

                return READER_SUCCESS;
	}
	return ln; // ln is READEDER ERROR as returned by this->findVariableLine, if not found
}
// end readInt
int Reader::readDlmLine(std::vector<double> &val, std::string dlm){
    val.clear();

    using std::cout;
    using std::endl;
    if (!file.good()) {
	cout << "bad file" << endl;
	return READER_BAD_FILE;
    }

    if (dlm.length()!= 1) {
	cout << "Error in Reader::readDlmLine - too many delimiters, fix me!\n" << endl;
	return READER_BAD_FORMAT;
    }
    std::string line;

    getline(file , line); // reads single line

    // REMOVE COMMENT LINES, WHITESPACE AT BEGINNING AND END OF FILE
    removeComments(line);
    removeBlanksAtBeginning(line);
    removeBlanksAtEnd(line);

    //size_t d = line.find_first_of(dlm);
    //cout << "delimiter at  " << d << endl;

    // EXTRAXT NUMBERS FROM LINE
    std::stringstream ss;
    ss << line;
    double a;
    while(ss >> a){
	val.push_back(a);
	//cout << a << endl;
    }

    return READER_SUCCESS;

}

int Reader::readNextAssignment(std::string& var , std::string& val){
	if ( ! this->file.good() )
		return READER_BAD_FILE;

	std::string line;	// READ SINGLE LINE FROM PREVIOUSLY OPENED FILE
	getline(file, line);

	var.clear();
	val.clear();

	this->removeComments(line);
	var = this->extractLHS(line);
	val = this->extractRHS(line);

	if ( (var.length() > 0) ) // IF var IS FOUND RETURN SUCCESS (val MAY BE EMPTY)
	    return READER_SUCCESS;

	return READER_NOT_FOUND;
}


int Reader::writeLine(const std::string& line, const int& ln){

	if ( ! file.good() )
		return READER_BAD_FILE;
	if (ln >= 0)
		if ( !this->gotoLine(ln) )
			file.clear();

	if ( file << line << endl )
		return READER_SUCCESS;
	else
		return READER_BAD_FILE;
}


bool Reader::removeBlanks(std::string& str){
    using std::string;
    string::iterator itr;
    itr = str.begin();

    string temp_str = "";
	string space    = " ";	// POSSIBLE "BLANK" CHARACTERS
    string tab      = "\t";
	string nl		= "\n";
	string ret		= "\r";
    for (size_t pos1 = 0 ; pos1 < str.length() ; pos1++){
		if ( (str.compare(pos1,1,space) != 0) &&
			 (str.compare(pos1,1,tab) != 0 ) &&
			 (str.compare(pos1,1,nl) != 0 ) &&
			 (str.compare(pos1,1,ret) != 0)	)
			 {
            temp_str.append( &str.at(pos1) , 1 ) ; // add 1 character only
            
        }
    }

    str.clear();
    str = temp_str;
return true;
}
//end removeBlanks
bool Reader::removeBlanksAtBeginning(std::string &str){
    using std::string;
    string::iterator itr;
    itr = str.begin();

    string space    = " ";
    string tab	    = "\t";
    size_t pos1;
    for (pos1 = 0 ; pos1 < str.length() ; pos1++){
	if ( (str.compare(pos1,1,space) == 0) || // compare characters 1-by-1
	     (str.compare(pos1,1,tab)) == 0) {
	  // if blank, do nothing
	}
	else
	{break;} // when non-blank found ->exit loop
    } // end for
    // copy non-blanks to a temp
    string temp;
    temp = str.substr(pos1 , str.length() - pos1);

    // replace original with temp
    str = temp;
    return true;
}

bool Reader::removeBlanksAtEnd(std::string &str){
	std::string whitespaces (" \t\f\v\n\r");
	size_t found;

	found=str.find_last_not_of(whitespaces);

	if (found!=std::string::npos)
	    str.erase(found+1);
	else
	    str.clear();            // str is all whitespace


    return true;

}

bool Reader::removeComments(std::string& str){

    size_t pos = str.find_first_of("#");
    str = str.substr(0 , pos);
    return true;
}
// end removeComments
void Reader::stringToLowercase(std::string& str){
    std::locale loc;
    std::string temp;
    temp = "";
    char c;
    for (unsigned int i = 0; i < str.length() ; i++){
	c = str[i];
	if (std::isupper(c , loc) )
	    c = std::tolower(c , loc);

	temp.append(1 , c);
	//std::cout << c << std::endl;
    }
    str = temp;

}
bool Reader::gotoLine(const unsigned int& i){

    file.seekg(0);
    std::string temp;
    for (unsigned int p = 0 ; p < i ; p++){
	getline(file, temp );
	if (file.eof())
	    return false;
    }
    return true;
}

std::string Reader::extractRHS( std::string line ){


	size_t eq = line.find('='); // location of '='
	if ( eq != std::string::npos ){ // if equals sign found
	    line = line.substr(eq+1 , std::string::npos); // keep everthing after '='
	    // cleanup
		removeComments(line);
	    removeBlanksAtBeginning(line);
	    removeBlanksAtEnd(line);
	    
		return line;
	}
        return std::string(""); // otherwise return empty string
}

std::string Reader::extractLHS(std::string line){
    size_t eq = line.find('=');
    if ( eq != std::string::npos){
	line = line.substr(0 , eq); // keep everything before '='
	removeBlanksAtBeginning(line);
	removeBlanksAtEnd(line);
	return line;
    }
    return std::string("");
}


int  Reader::stringToDoubleVec(std::string& str , std::vector <double>& vec){
	vec.clear();
	size_t opn = str.find('[');
	size_t cls = str.find(']');
	// ERROR CHECKING. MAY NEED MORE
	if ( (cls < opn) || 				// if closing bracket first
		 (opn == std::string::npos) ||	// if open bracket not found
		 (cls == std::string::npos) )   // if closing bracket not found
		 	return READER_BAD_FORMAT;

	str = str.substr(opn+1, cls-1); // REMOVE BRACKETS


	// START PUSHING COMMA SEPARATED DOUBLES ONTO A VECTOR
	size_t pos = str.find_first_of(',');
	double temp = 0;

	while ( (pos < std::string::npos ) ){

		if (! (std::stringstream( str.substr(0,pos) ) >> temp ) ) // STRINGSTREAM TO TEST CHARACTER VALIDITY
			return READER_BAD_FORMAT;

		vec.push_back(temp);
		str = str.substr(pos+1 , std::string::npos );
		pos = str.find_first_of(',');
	}

	if (! (std::stringstream( str.substr(0,pos) ) >> temp ) )	// LAST ITEM IN LIST NEEDS TO DO SEPARATELY
			return READER_BAD_FORMAT;
	vec.push_back(temp);

	return READER_SUCCESS;

}



#include <reader.h>

Reader::Reader(){
	//file = new std::fstream();

}

Reader::~Reader(){
	//delete file;
}

bool Reader::openFile(const std::string& filename){

	//file.open( filename.c_str() );
	file.open( filename.c_str() );
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
      //  case READER_BAD_FILE	:{
//	    return string("READER_BAD_FILE");
//	    break;
//	}
        default:{
            return string("UNDEFINED ERROR");
            break;
            }
        }

}
//end getErrorString

int Reader::findStringLine(std::string str,
			   unsigned int offset
			   ){
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


int Reader::readString(std::string var, std::string& val){
    // std::cout << "searching for :" << var << std::endl;
	file.seekg(0);
	if (file.good() ) {

		std::string line;

	while( ! file.eof() ){
		getline(file , line);

        removeBlanks(line);     // remove extra rubbish
        removeComments(line);

        if ( line.find(var) == 0 )
        {
            size_t pos = line.find ("="); // find "=" character
            if (pos == std::string::npos){
                std::cout << "line is:" << line << std::endl;
                return READER_BAD_FORMAT;
            }

            std::string s_val = line.substr(pos+1 , line.length() - (pos+1) );

            //removeBlanks(s_val);
            val = s_val;
            return READER_SUCCESS;
        }
    }
	file.seekg(0);
	file.clear();  // clear eof bit;
    return READER_NOT_FOUND;
    }
    std::cout << "bad file" << std::endl;
    return READER_NOT_FOUND;//READER_BAD_FILE;

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


bool Reader::removeBlanks(std::string& str){
    using std::string;
    string::iterator itr;
    itr = str.begin();

    string temp_str = "";
    string space    = " ";
    string tab      = "\t";
    for (size_t pos1 = 0 ; pos1 < str.length() ; pos1++){
        if ( (str.compare(pos1,1,space) != 0) && ( str.compare(pos1,1,tab) != 0 ) ){
            temp_str.append( &str.at(pos1) , 1 ) ; // add 1 character only
            //std::cout << "character is " << str.at(pos1) << std::endl;
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

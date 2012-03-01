#include <filesysfun.h>

namespace FilesysFun
{

bool setCurrentDirectory(const std::string& destdir)
{
// Changes current working directory to destdir.
// returns true/false if/not succesful.
#ifdef Linux
    if ( chdir( destdir.c_str() ) == 0 )
        return true;
    else
        return false;
#endif

#ifdef Windows

    bool success = SetCurrentDirectoryA( destdir.c_str() );

    return success;
#endif

}

std::string getCurrentDirectory()
{
    std::string cwd;

#ifdef Linux
    char* buff = get_current_dir_name(); // allocates space to buffer
    cwd = buff;
    free( buff );   // must free buffer
#endif

#ifdef Windows
    char buff[FILENAME_MAX];
    _getcwd( buff, sizeof(buff) );
    cwd = buff;
#endif

    return cwd;

}

bool dirExists(const std::string& dir)
{
    // checks whether directory dir exists, returning true/false

    std::string curdir = getCurrentDirectory();

    if ( setCurrentDirectory( dir ) ) // if can go to queried dir
    {
        setCurrentDirectory( curdir ); // go back to previous
        return true;
    }

    return false; // otherwise assume dir does not exist
}

bool fileExists(const std::string &file)
{
    // CHECK WHETHER FILE EXISTS. RETURN TRUE/FALSE

    // TRY TO OPEN FILE. IF THIS FAILS RETURN FALSE
    FILE* fid;
    fid = fopen(file.c_str(), "r" );

    if ( fid )
    {
        fclose(fid);
        return true;
    }

    return false;
}

bool createDirectory(const std::string& newdir)
{
#ifdef Linux
    // set directory permission, read, write, execute
    // for User, Group and Other
    mode_t perm =   S_IRWXU | S_IRWXG | S_IRWXO;

    if( mkdir( newdir.c_str(), perm ) != 0 )
        return false;
    else
        return true; //success
#endif

#ifdef Windows
    return CreateDirectoryA( newdir.c_str(), NULL );
#endif
}

bool copyFile(const std::string &srcFile,
              const std::string &dstDir,
              const std::string &dstFile)
{
    // COPIES SOURCE FILE TO DESTINATION FILE

    if ( !fileExists(srcFile) )
    {
        printf("error in %s %s, source file %s does not exist\n",__func__, __FILE__, srcFile.c_str() );
        return false;
    }

    std::ifstream srcF(srcFile.c_str(), std::fstream::in | std::fstream::binary);

    // MAKE SURE OPENEND OK
    if (!srcF.good() )
    {
        printf("error in %s, %s, could not open source file %s\n", __func__, __FILE__, srcFile.c_str() );
        srcF.close();
        return false;
    }

    // CHECK THAT OUTPUT DIR EXISTS
    if ( !dirExists(dstDir) )
    {
        if (!createDirectory(dstDir) )
        {
            printf("error in %s, %s, could not create destination directory %s\n", __func__,
                   __FILE__, dstDir.c_str() );
            srcF.close();
            return false;
        }
    }


    // OPEN OUTPUT FILE
    std::string fout = dstDir + "/";
    if (dstFile.length() == 0 )
        fout += srcFile;
    else
        fout += dstFile;

    std::ofstream dstF(fout.c_str(), std::fstream::out |std::fstream::binary);
    // CHECK OPENEED OK
    if ( !dstF.good() )
    {
        printf("error in %s, %s, could not create destination file %s\n",
               __func__, __FILE__, fout.c_str() );
        srcF.close();
        dstF.close();
        return false;
    }


    // WRITE DATA FROM SOURCE TO DESTINATION
    dstF << srcF.rdbuf();


    srcF.close();
    dstF.close();
    return true;
}


}//end namespace FilesysFun

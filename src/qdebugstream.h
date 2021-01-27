/* IO redirecting 
 * 1. http://www.qtforum.org/article/39768/redirecting-std-cout-std-cerf-qdebug-to-qtextedit.html
 * 2. http://stackoverflow.com/questions/12978973/unknown-ouput-with-qdebugstream-and-qtextedit
 * 3. http://stackoverflow.com/questions/10308425/redirect-stdcout-to-a-qtextedit
 * 4. http://stackoverflow.com/questions/9211298/redirect-stdcout-to-qtextedit?answertab=oldest#tab-top
 * 5. http://www.archivum.info/qt-interest@trolltech.com/2005-06/00040/Redirecting-std-cout-to-QTextEdit--Using-QTextEdit-as-a-Log-Window.html
 * Where 5 was half implemented then removed until later. 
 * Below code copied then modified from website 5. */


#ifndef QDEBUGSTREAM_H
#define	QDEBUGSTREAM_H

//################
//# qdebugstream.h  #
//################

#include <iostream>
#include <streambuf>
#include <string>

//#include "qtextedit.h"

class QPlainTextEdit;

class QDebugStream : public std::basic_streambuf<char>
{
public:
 QDebugStream(std::ostream &stream, QPlainTextEdit* text_edit);// : m_stream(stream);
 ~QDebugStream();

protected:
 virtual int_type overflow(int_type v);
 virtual std::streamsize xsputn(const char *p, std::streamsize n);

private:
 std::ostream &m_stream;
 std::streambuf *m_old_buf;
 std::string m_string;
 QPlainTextEdit* log_window;
};

#endif	/* QDEBUGSTREAM_H */


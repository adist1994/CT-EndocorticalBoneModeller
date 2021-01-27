#include <QtWidgets>
#include <iostream>
#include <streambuf>
#include <string>

#include "qdebugstream.h"

QDebugStream::QDebugStream(std::ostream &stream, QPlainTextEdit* text_edit) : m_stream(stream)
 {
  log_window = text_edit;
  m_old_buf = stream.rdbuf();
  stream.rdbuf(this);
 }

 QDebugStream::~QDebugStream()
 {
  // output anything that is left
  if (!m_string.empty()) {
    //log_window->appendPlainText(m_string.c_str());
    log_window->appendHtml(m_string.c_str());
  }

  m_stream.rdbuf(m_old_buf);
 }

 std::streambuf::int_type QDebugStream::overflow(int_type v)
 {
  if (v == '\n')
  {
   //log_window->appendPlainText(m_string.c_str());
   log_window->appendHtml(m_string.c_str());
   m_string.erase(m_string.begin(), m_string.end());
  }
  else
   m_string += v;

  return v;
 }

 std::streamsize QDebugStream::xsputn(const char *p, std::streamsize n)
 {
  m_string.append(p, p + n);

  int pos = 0;
  while (pos != std::string::npos)
  {
   pos = m_string.find('\n');
   if (pos != std::string::npos)
   {
    std::string tmp(m_string.begin(), m_string.begin() + pos);
    //log_window->appendPlainText(tmp.c_str());
    log_window->appendHtml(m_string.c_str());
    m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
   }
  }

  return n;
 }

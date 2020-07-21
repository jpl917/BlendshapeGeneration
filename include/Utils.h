#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <string>
using namespace std;

#if DEBUG_BUILD
#define DEBUG_OUTPUT(x) { cout << x << endl; }
#define DEBUG_EXPR(x) {x}
#else
#define DEBUG_OUTPUT(x)
#define DEBUG_EXPR(x)
#endif

namespace ColorOutput {
  static const string Black = "\033[0;30m";
  static const string Red = "\033[0;31m";
  static const string Green = "\033[0;32m";
  static const string Orange = "\033[0;33m";
  static const string Blue = "\033[0;34m";
  static const string Purple = "\033[0;35m";
  static const string Cyan = "\033[0;36m";
  static const string LightGray = "\033[0;37m";

  static const string DarkBlack = "\033[1;30m";
  static const string LightRed = "\033[1;31m";
  static const string LightGreen = "\033[1;32m";
  static const string Yellow = "\033[1;33m";
  static const string LightBlue = "\033[1;34m";
  static const string LightPurple = "\033[1;35m";
  static const string LightCyan = "\033[1;36m";
  static const string White = "\033[1;37m";

  static const string Reset = "\033[0m";
};

class ColorStream {
public:
  ColorStream(const string& color, ostream& out = cout) : out(out){
    out << color;
  }
  ~ColorStream() {
    out << ColorOutput::Reset << endl;
  }
  template <typename T>
  ColorStream& operator<<(const T& t) {
    out << t;
    return (*this);
  }

private:
  ostream& out;
};


// debugging related
// dummy
static void debug(){}

template <typename T>
void debug(const string& name, T value)
{
  cout << name << " = " << value << endl;
}

template <typename T, typename ...Args>
void debug(const string& name, T value, Args ...args)
{
  cout << name << " = " << value << endl;
  debug(args...);
}

// general console message
inline void message(const string& msg) {
  cout << msg << endl;
}

inline void error(const string& msg) {
  cerr << "Error:\t" << msg << endl;
}

inline void abort(const string& msg) {
  cerr << "Critical error:\t" << msg << endl;
  exit(0);
}

#endif

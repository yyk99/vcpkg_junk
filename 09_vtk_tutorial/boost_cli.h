//
//
//

#ifndef BOOST_CLI_H
#define BOOST_CLI_H

#include <boost/program_options.hpp>
#include <string>

namespace po = boost::program_options;

namespace CLI{
class App {
public:
    App() {}
    void add_flag("-w, --wireframe", wireframeOn, "Render a wireframe.");
};
}
#endif
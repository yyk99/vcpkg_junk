//
//
//

#include <proj.h>

#include <gtest/gtest.h>

#include "../common/DebuggingConsole.h"

class ProjF : public testing::Test {

};

TEST_F(ProjF, version) { CONSOLE("Hello..."); }

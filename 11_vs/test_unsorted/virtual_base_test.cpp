//
//
//

#include <gtest/gtest.h>
#include <stdint.h>

#include <cstddef>
#include <string>

#include "test_unsorted.h"

class Widget {
    void *handle;

public:
    Widget() : handle(0) {}
};

class Data {
public:
    int32_t property1;

public:
    Data() : property1{} {}

    virtual ~Data() {}
    virtual const char *type_name() const { return "class Data"; }
};

#include "virtual_base_test_plain.h"
#include "virtual_base_test_virt.h"

class VirtualBaseF : public testing::Test {};

TEST_F(VirtualBaseF, t1) {
    using namespace plain;

    CONSOLE_EVAL(sizeof(Widget));
    CONSOLE_EVAL(sizeof(Data));

    CONSOLE_EVAL(sizeof(Dialog));
    CONSOLE_EVAL(sizeof(DataMore));
    CONSOLE_EVAL(sizeof(DialogMore));

    Dialog d1;
    EXPECT_EQ(0, d1.property1);
    EXPECT_STREQ("class Data", d1.type_name());

    DialogMore d2;
    EXPECT_STREQ("class Data", d2.Dialog::type_name());
#if _MSC_VER && !defined(__clang__)
    EXPECT_STREQ("class Data", d2.Data::type_name());
#endif
    EXPECT_STREQ("class DataMore", d2.DataMore::type_name());

#if _MSC_VER && !defined(__clang__)
    EXPECT_EQ(0, d2.Data::property1);
#endif
    EXPECT_EQ(0, d2.Dialog::property1);
    EXPECT_EQ(42, d2.DataMore::property1);
    EXPECT_EQ(42, d2.DataMore::property2);

    Process p;
    p.process(d1);
    p.process(static_cast<DataMore &>(d2));
    p.process(static_cast<Dialog &>(d2));
}

TEST_F(VirtualBaseF, t2) {
    using namespace virt;

    CONSOLE_EVAL(sizeof(Widget));
    CONSOLE_EVAL(sizeof(Data));

    CONSOLE_EVAL(sizeof(Dialog));
    CONSOLE_EVAL(sizeof(DataMore));
    CONSOLE_EVAL(sizeof(DialogMore));

    Dialog d1;
    EXPECT_EQ(0, d1.property1);
    EXPECT_STREQ("class Data", d1.type_name());

    DialogMore d2;
    EXPECT_EQ(42, d2.Data::property1);
    EXPECT_EQ(42, d2.DataMore::property1);
    EXPECT_EQ(42, d2.DataMore::property2);

    EXPECT_EQ(42, d2.property1);
    EXPECT_EQ(42, d2.property2);

    Process p;
    p.process(d1);
    p.process(d2);
    p.process(static_cast<DataMore &>(d2));
    p.process(static_cast<Dialog &>(d2));
}

struct base {
    virtual int number() { return 0; }
};
struct weak : public virtual base {
    void print() { // seems to only depend on base, but depends on dominant
        std::cout << number() << std::endl;
    }
};
struct dominant : public virtual base {
    int number() override { return 5; }
};
struct derived : public weak, public dominant {};

TEST_F(VirtualBaseF, t3) {
    weak w;
    w.print(); // 0
    derived d;
    d.print(); // 5
}

/* One more design. The data is stored in pointers */

namespace ptr {

class Dialog : public Widget {
public:
    std::unique_ptr<Data> m_data;
    virtual Data &data() { return *m_data.get(); };

public:
    Dialog() : m_data{new Data}, control1(120) {}

    int32_t control1;
};

class DataMore : public Data {
public:
    int32_t property2;

public:
    DataMore() : property2{42} { property1 = 42; }

    const char *type_name() const override { return "class DataMore"; }
};

class DialogMore : public Dialog {
public:
    DialogMore() : control2(55) { m_data.reset(new DataMore); }

    int32_t control2;

    DataMore &data() override {
        return *dynamic_cast<DataMore *>(m_data.get());
    };
};

class Process {
public:
    void process(Data const &data) {
        std::cout << "data type: " << data.type_name() << ":" << data.property1
                  << std::endl;
    }
};

} // namespace ptr

TEST_F(VirtualBaseF, t4) {
    using namespace ptr;

    CONSOLE_EVAL(sizeof(Widget));
    CONSOLE_EVAL(sizeof(Data));

    CONSOLE_EVAL(sizeof(Dialog));
    CONSOLE_EVAL(sizeof(DataMore));
    CONSOLE_EVAL(sizeof(DialogMore));

    Dialog d1;
    EXPECT_EQ(0, d1.data().property1);
    EXPECT_STREQ("class Data", d1.data().type_name());

    DialogMore d2;
    EXPECT_EQ(42, d2.data().property1);
    EXPECT_EQ(42, d2.data().property1);
    EXPECT_EQ(42, d2.data().property2);

    EXPECT_EQ(42, d2.data().property1);
    EXPECT_EQ(42, d2.data().property2);

    Dialog &dp1 = d2;

    dp1.data().property1 = 33;
    EXPECT_EQ(33, d2.data().property1);

    Process p;
    p.process(d1.data());
    p.process(d2.data());
}

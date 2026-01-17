#pragma once

namespace virt {

class Dialog : public Widget, public virtual Data {
public:
    Dialog() : control1(120) {}

    int32_t control1;
};

class DataMore : public virtual Data {
public:
    int32_t property2;

public:
    DataMore() : property2{42} { property1 = 42; }

    const char *type_name() const override { return "class DataMore"; }
};

class DialogMore : public Dialog, public DataMore {
public:
    DialogMore() : control2(55) {}

    int32_t control2;
};

class Process {
public:
    void process(Data const &data) {
        std::cout << "data type: " << data.type_name() << ":" << data.property1
                  << std::endl;
    }
};

} // namespace virt

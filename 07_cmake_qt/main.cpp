//
//
//
#include <QApplication>
#include <QMessageBox>
#include <QPushButton>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QPushButton button("Hello, Qt6!");

    QAction aboutAction("About", &button);
    QObject::connect(&aboutAction, &QAction::triggered,
                     [&]() { QMessageBox::aboutQt(&button, "About"); });
    button.addAction(&aboutAction);
    button.setContextMenuPolicy(Qt::ActionsContextMenu);
    button.resize(200, 60);
    button.show();

    return app.exec();
}
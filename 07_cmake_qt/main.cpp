//
//
//
#include <QApplication>
#include <QPushButton>
#include <QMessageBox>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QPushButton button("Hello, Qt6!");

    QAction aboutAction("About", &button);
    QObject::connect(&aboutAction, &QAction::triggered, [&]() {
        QMessageBox::about(&button, "About", "This is a Qt6 application.");
    });
    button.addAction(&aboutAction);
    button.setContextMenuPolicy(Qt::ActionsContextMenu);
    button.resize(200, 60);
    button.show();

    return app.exec();
}
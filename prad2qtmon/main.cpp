// ─────────────────────────────────────────────────────────────────────────────
// prad2_monitor – Lightweight Qt WebEngine client for PRad2 event viewer/monitor
//
// Usage:
//   prad2_monitor                        # http://localhost:5051
//   prad2_monitor -H clonpc19            # http://clonpc19:5051
//   prad2_monitor -H clonpc19 -p 8080   # http://clonpc19:8080
// ─────────────────────────────────────────────────────────────────────────────

#include <QApplication>
#include <QWebEngineView>
#include <QWebEnginePage>
#include <QUrl>

#include <iostream>
#include <string>
#include <unistd.h>

static void printUsage(const char *prog)
{
    std::cerr
        << "Usage:\n"
        << "  " << prog << " [-H host] [-p port]\n"
        << "\nOptions:\n"
        << "  -H <host>   Server hostname (default: localhost)\n"
        << "  -p <port>   Server port (default: 5051)\n"
        << "  -h          Show this help\n";
}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("PRad2 Monitor");
    app.setApplicationVersion("1.0.0");

    std::string host = "localhost";
    std::string port = "5051";

    optind = 1;
    int opt;
    while ((opt = getopt(argc, argv, "H:p:h")) != -1) {
        switch (opt) {
        case 'H': host = optarg; break;
        case 'p': port = optarg; break;
        case 'h': printUsage(argv[0]); return 0;
        default:  printUsage(argv[0]); return 1;
        }
    }

    QUrl url(QString("http://%1:%2")
        .arg(QString::fromStdString(host), QString::fromStdString(port)));

    std::cout << "Loading: " << url.toString().toStdString() << "\n";

    QWebEngineView view;
    view.setWindowTitle("PRad2 Monitor — " + url.toString());
    view.resize(1280, 800);
    // log load errors to console
    QObject::connect(view.page(), &QWebEnginePage::loadFinished,
                     [&url](bool ok) {
        if (!ok)
            std::cerr << "Failed to load: " << url.toString().toStdString() << "\n";
    });

    view.setUrl(url);
    view.show();

    return app.exec();
}

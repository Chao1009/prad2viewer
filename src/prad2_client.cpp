// ─────────────────────────────────────────────────────────────────────────────
// prad2_client – Lightweight Qt WebEngine client for a remote PRad2 server
//
// Connects to a running prad2_server instance (file viewer or online monitor).
//
// Usage:
//   prad2_client                        # http://localhost:5051
//   prad2_client -H clonpc19            # http://clonpc19:5051
//   prad2_client -H clonpc19 -p 8080   # http://clonpc19:8080
// ─────────────────────────────────────────────────────────────────────────────

#include <QApplication>
#include <QWebEngineView>
#include <QWebEnginePage>
#include <QUrl>

#include <iostream>
#include <string>
#include <getopt.h>

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
    app.setApplicationName("PRad2 Client");
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

    std::cout << "Connecting to: " << url.toString().toStdString() << "\n";

    QWebEngineView view;
    view.setWindowTitle("PRad2 Client \u2014 " + url.toString());
    view.resize(1280, 800);
    QObject::connect(view.page(), &QWebEnginePage::loadFinished,
                     [&url](bool ok) {
        if (!ok)
            std::cerr << "Failed to load: " << url.toString().toStdString() << "\n";
    });

    view.setUrl(url);
    view.show();

    return app.exec();
}

// =========================================================================
// prad2evviewer — Standalone Qt viewer for PRad-II events
//
// Embeds a ViewerServer on localhost and displays the web frontend via
// Qt WebEngine. All updates to the web frontend or server API are
// automatically reflected here — zero maintenance.
//
// Usage:
//   prad2evviewer                          # empty, use File > Open
//   prad2evviewer data.evio                # open file directly
//   prad2evviewer data.evio -H             # open with histograms
//   prad2evviewer -d /data/stage6          # enable file browser
// =========================================================================

#include "viewer_server.h"

#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>
#include <QStatusBar>
#include <QFileDialog>
#include <QFileInfo>
#include <QAction>
#include <QWebEngineView>
#include <QWebEnginePage>
#include <QTimer>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QUrl>
#include <QLabel>

#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif
#ifndef RESOURCE_DIR
#define RESOURCE_DIR "."
#endif

// ─────────────────────────────────────────────────────────────────────────
// Main window
// ─────────────────────────────────────────────────────────────────────────

class ViewerWindow : public QMainWindow {
    Q_OBJECT
public:
    ViewerWindow(ViewerServer *server, bool histDefault)
        : server_(server), histDefault_(histDefault)
    {
        setWindowTitle("PRad-II Event Viewer");
        resize(1360, 860);
        setAcceptDrops(true);

        // --- WebEngine view ---
        view_ = new QWebEngineView(this);
        setCentralWidget(view_);

        QUrl url(QString("http://localhost:%1").arg(server_->port()));
        view_->setUrl(url);

        // --- Menu bar ---
        auto *fileMenu = menuBar()->addMenu("&File");
        auto *openAct = fileMenu->addAction("&Open...", this, &ViewerWindow::onOpen,
                                            QKeySequence::Open);
        Q_UNUSED(openAct);
        fileMenu->addAction("Open with &Histograms...", this,
                            &ViewerWindow::onOpenHist);
        fileMenu->addSeparator();
        fileMenu->addAction("E&xit", this, &QWidget::close, QKeySequence::Quit);

        auto *viewMenu = menuBar()->addMenu("&View");
        viewMenu->addAction("&Reload", this, [this]() { view_->reload(); },
                            QKeySequence::Refresh);
#ifdef WITH_ET
        etAction_ = viewMenu->addAction("Go &Online", this,
                                        &ViewerWindow::onToggleOnline);
#endif

        // --- Status bar ---
        statusLabel_ = new QLabel("Ready");
        statusBar()->addWidget(statusLabel_, 1);

        // --- Progress timer ---
        progressTimer_ = new QTimer(this);
        connect(progressTimer_, &QTimer::timeout, this,
                &ViewerWindow::pollProgress);
        progressTimer_->start(500);
    }

protected:
    void dragEnterEvent(QDragEnterEvent *e) override {
        if (e->mimeData()->hasUrls()) e->acceptProposedAction();
    }
    void dropEvent(QDropEvent *e) override {
        for (auto &url : e->mimeData()->urls()) {
            QString path = url.toLocalFile();
            if (path.contains(".evio")) {
                openFile(path.toStdString(), histDefault_);
                break;
            }
        }
    }

private slots:
    void onOpen() {
        QString path = QFileDialog::getOpenFileName(
            this, "Open EVIO File", lastDir_,
            "EVIO Files (*.evio *.evio.*);;All Files (*)");
        if (!path.isEmpty()) {
            lastDir_ = QFileInfo(path).absolutePath();
            openFile(path.toStdString(), false);
        }
    }
    void onOpenHist() {
        QString path = QFileDialog::getOpenFileName(
            this, "Open EVIO File (with histograms)", lastDir_,
            "EVIO Files (*.evio *.evio.*);;All Files (*)");
        if (!path.isEmpty()) {
            lastDir_ = QFileInfo(path).absolutePath();
            openFile(path.toStdString(), true);
        }
    }
#ifdef WITH_ET
    void onToggleOnline() {
        if (server_->mode() == "online") {
            view_->page()->runJavaScript(
                "fetch('/api/mode/file',{method:'POST'})");
        } else {
            // open the ET connection dialog in the web frontend
            view_->page()->runJavaScript("openEtDialog()");
        }
    }
#endif
    void pollProgress() {
        if (server_->isLoading()) {
            auto p = server_->getProgress();
            QString phase = QString::fromStdString(
                p.value("phase", "loading"));
            int cur = p.value("current", 0);
            int tot = p.value("total", 0);
            statusLabel_->setText(
                QString("Loading: %1 (%2/%3)")
                    .arg(phase).arg(cur).arg(tot));
        } else {
            auto m = QString::fromStdString(server_->mode());
            statusLabel_->setText(QString("Mode: %1").arg(m));
        }
#ifdef WITH_ET
        if (etAction_) {
            etAction_->setText(server_->mode() == "online"
                               ? "Go to &File Viewer"
                               : "Go &Online");
        }
#endif
    }

private:
    void openFile(const std::string &path, bool hist) {
        server_->loadFile(path, hist);
        setWindowTitle(QString("PRad-II Event Viewer \u2014 %1")
                           .arg(QFileInfo(QString::fromStdString(path)).fileName()));
        // reload page so frontend picks up the new file
        QTimer::singleShot(200, this, [this]() { view_->reload(); });
    }

    ViewerServer *server_;
    QWebEngineView *view_;
    QLabel *statusLabel_;
    QTimer *progressTimer_;
    bool histDefault_;
    QString lastDir_;
#ifdef WITH_ET
    QAction *etAction_ = nullptr;
#endif
};

// ─────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("PRad-II Event Viewer");
    app.setApplicationVersion("1.0.0");

    ViewerServer::Config cfg;
    cfg.database_dir = DATABASE_DIR;
    cfg.resource_dir = RESOURCE_DIR;
    cfg.port = 0; // auto-select port
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  cfg.database_dir = env;
    if (const char *env = std::getenv("PRAD2_RESOURCE_DIR"))  cfg.resource_dir = env;

    static struct option long_opts[] = {
        {"hist",       no_argument,       nullptr, 'H'},
        {"data-dir",   required_argument, nullptr, 'd'},
        {"config",     required_argument, nullptr, 'c'},
        {"daq-config", required_argument, nullptr, 'D'},
        {"help",       no_argument,       nullptr, '?'},
        {nullptr, 0, nullptr, 0},
    };

    optind = 1;
    int opt;
    while ((opt = getopt_long(argc, argv, "Hd:c:D:", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'H': cfg.hist_enabled = true; break;
        case 'd': cfg.data_dir = optarg; break;
        case 'c': cfg.config_file = optarg; break;
        case 'D': cfg.daq_config_file = optarg; break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " [evio_file] [-H] [-d data_dir]"
                      << " [-c config.json] [-D daq_config.json]\n";
            return 1;
        }
    }
    if (optind < argc) cfg.initial_file = argv[optind];

    // init & start server
    ViewerServer server;
    server.init(cfg);
    int port = server.startAsync();

    std::cout << "Viewer at http://localhost:" << port << "\n";

    // create window
    ViewerWindow window(&server, cfg.hist_enabled);
    window.show();

    int rc = app.exec();

    server.stop();
    return rc;
}

#include "prad2evviewer.moc"

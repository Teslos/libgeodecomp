#ifndef _libgeodecomp_examples_latticegas_framegrabber_h_
#define _libgeodecomp_examples_latticegas_framegrabber_h_

#include <iostream>
#include <stdexcept>
#include <QObject>
#include <libgeodecomp/examples/latticegas/fpscounter.h>

class FrameGrabber : public QObject, FPSCounter
{
    Q_OBJECT

public:
    FrameGrabber(QObject *parent);

    ~FrameGrabber();

public slots:
    void grab();

    void info()
    {
        std::cout << "FrameGrabber @ " << fps() << " FPS\n";
    }

signals:
    void newFrame(unsigned *frame);

private:
    // ugly hack. we can't include cv.h here since nvcc won't compile it
    void *capture;
};

#endif

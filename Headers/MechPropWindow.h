#ifdef _USE_QT_

#ifndef MECHPROPWINDOW_H
#define MECHPROPWINDOW_H

#include <QWidget>
#include <QDialog>
#include <QSpinBox>
#include <QGraphicsView>
#include "SpecialQSpinBox.h"
#include "SpecialQDoubleSpinBox.h"
#include "Tissue.h"

#define PI 3.14159265

QT_BEGIN_NAMESPACE
class QGroupBox;
class QLabel;
class QLineEdit;
class MainWindow;
QT_END_NAMESPACE


class MechPropWindow : public QWidget
{
    Q_OBJECT
public:
    MechPropWindow(MainWindow * mainWindow);
    
signals:
    //void signal_generateTissue( int kind, bool isPeriodic=false, int noRows=0, int noCols=0, double height=0, double width=0, double length=0 );

private:
    //MechanicalTissueProperties mechProp;
    
    
private slots:
    void slot_show();
    //void slot_Ta_1(double x){mechProp.3}
};


#include "MainWindow.h"

#endif

#endif // _USE_QT_
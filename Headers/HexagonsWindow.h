#ifdef _USE_QT_

#ifndef HexagonsWindow_H
#define HexagonsWindow_H

#include <QWidget>
#include <QDialog>
#include <QSpinBox>
#include <QCheckBox>
#include <QGraphicsView>
#include <QtGui>
#include "SpecialQSpinBox.h"
#include "SpecialQDoubleSpinBox.h"

QT_BEGIN_NAMESPACE
class QGroupBox;
class QLabel;
class QLineEdit;
class MainWindow;
QT_END_NAMESPACE


class HexagonsWindow : public QWidget
{
    Q_OBJECT
public:
    HexagonsWindow(MainWindow * mainWindow);
    
private:
    bool isPeriodic;
    int noRows;
    int noCols;
    double height;
    double width;
    double length;
    QCheckBox *QC_isPeriodic;
    QPushButton *QP_create_round_tissue, *QP_generate;
    SpecialQSpinBox *QS_noCols, *QS_noRows, *QS_radiusInCells;
    SpecialQDoubleSpinBox *QS_height, *QS_width, *QS_length;
    
signals:
    void signal_generateTissue( int kind, bool isPeriodic=false, int noRows=0, int noCols=0, double height=0, double width=0, double length=0 );
    void signal_createRoundHexagonalTissue(int, double, double, double);

private slots:
    void slot_isPeriodic(int state){isPeriodic=(state==2);}
    void slot_height(double h){height=h;}
    void slot_width(double w){width=w;}
    void slot_length(double l){length=l;}
    void slot_noRows(int nr){noRows=nr;}
    void slot_noCols(int nc){noCols=nc;}
    void slot_generate();
    void slot_generate_cyst();
    void slot_generate_bigCyst();
    void slot_show();
    void slot_createRoundHexagonalTissue();
};


#include "MainWindow.h"

#endif
#endif // _USE_QT_
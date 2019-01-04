#ifdef _USE_QT_

#ifndef SpecialQDoubleSpinBox_H
#define SpecialQDoubleSpinBox_H

#include <QApplication>
#include <QWidget>
#include <QGraphicsView>
#include<QString>
#include<QSpinBox>
#include<QEvent>


class SpecialQDoubleSpinBox: public QDoubleSpinBox
{
    Q_OBJECT
    public :
    SpecialQDoubleSpinBox(QWidget* parent=0);
    //bool eventFilter(QObject *anObject, QEvent *anEvent);
    
signals:
    void ValueChangedAfterEnter(double);
    public slots:
    void EditingFinished();
    
};
#endif

#endif // _USE_QT_
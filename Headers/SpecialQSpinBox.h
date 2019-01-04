#ifdef _USE_QT_

#ifndef SpecialQSpinBox_H
#define SpecialQSpinBox_H

#include <QApplication>
#include <QWidget>
#include <QGraphicsView>
#include <QString>
#include <QSpinBox>
#include <QEvent>


class SpecialQSpinBox: public QSpinBox
{
	Q_OBJECT
	public :
	SpecialQSpinBox(QWidget* parent=0);
    //bool eventFilter(QObject *anObject, QEvent *anEvent);
	
	signals:
	void ValueChangedAfterEnter(int);
public slots:
	void EditingFinished();
};
#endif

#endif // _USE_QT_
#ifdef _USE_QT_


#include "SpecialQSpinBox.h"


SpecialQSpinBox::SpecialQSpinBox(QWidget *parent)
: QSpinBox(parent)
{
	//QinstallEventFilter();
	connect(this,SIGNAL(editingFinished()),this,SLOT(EditingFinished()));
};

void SpecialQSpinBox::EditingFinished(){
emit(ValueChangedAfterEnter(value()));
}
/*bool SpecialQSpinBox::eventFilter(QObject *anObject, QEvent *anEvent)
{
	switch (anEvent->type())
    {
		case QEvent::KeyPress:
		{
			if((anEvent->key() == Qt::Key_Return) || (anEvent->key() == Qt::Key_Enter))
			{
				emit ValueChangedAfterEnter(value());
			}
			break;
		}
    }
	return QSpinBox::eventFilter(anObject, anEvent);
}*/

#endif // _USE_QT_
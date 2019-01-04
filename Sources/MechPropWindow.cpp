#ifdef _USE_QT_

#include <QtGui>
#include "MechPropWindow.h"

MechPropWindow::MechPropWindow(MainWindow * mainWindow)
: QWidget()
{
    /* DEFINE THE GEOMETRY */
	//setGeometry(40, 40, 500, 500);
    setWindowTitle("Mechanical properties");

    /* CREATE ALL ELEMENTS*/
    // headlines of the tables
    QFont QF_headline("Times", 12, QFont::Bold);
    QF_headline.setUnderline(true);
    QLabel *QL_T_l = new QLabel("T_a");
    QL_T_l->setFont(QF_headline);
    QLabel *QL_G_a = new QLabel("G_a");
    QL_G_a->setFont(QF_headline);
    QLabel *QL_G_b = new QLabel("G_b");
    QL_G_b->setFont(QF_headline);
    
    // labels
    QFont QF_normal("Times", 10);
    // first matrix: lateral surface tensions
    QLabel *QL_noCell_1 = new QLabel("No Cell");
    QL_noCell_1->setFont(QF_normal);
    QLabel *QL_type1_v_1 = new QLabel("Type 1");
    QL_type1_v_1->setFont(QF_normal);
    QLabel *QL_type2_v_1 = new QLabel("Type 2");
    QL_type2_v_1->setFont(QF_normal);
    QLabel *QL_type3_v_1 = new QLabel("Type 3");
    QL_type3_v_1->setFont(QF_normal);
    QLabel *QL_type1_h_1 = new QLabel("Type 1");
    QL_type1_h_1->setFont(QF_normal);
    QLabel *QL_type2_h_1 = new QLabel("Type 2");
    QL_type2_h_1->setFont(QF_normal);
    QLabel *QL_type3_h_1 = new QLabel("Type 3");
    QL_type3_h_1->setFont(QF_normal);
    
    // second matrix: apical line tensions
    QLabel *QL_noCell_2 = new QLabel("No Cell");
    QL_noCell_2->setFont(QF_normal);
    QLabel *QL_type1_v_2 = new QLabel("Type 1");
    QL_type1_v_2->setFont(QF_normal);
    QLabel *QL_type2_v_2 = new QLabel("Type 2");
    QL_type2_v_2->setFont(QF_normal);
    QLabel *QL_type3_v_2 = new QLabel("Type 3");
    QL_type3_v_2->setFont(QF_normal);
    QLabel *QL_type1_h_2 = new QLabel("Type 1");
    QL_type1_h_2->setFont(QF_normal);
    QLabel *QL_type2_h_2 = new QLabel("Type 2");
    QL_type2_h_2->setFont(QF_normal);
    QLabel *QL_type3_h_2 = new QLabel("Type 3");
    QL_type3_h_2->setFont(QF_normal);
    
    // third matrix: basal line tensions
    QLabel *QL_noCell_3 = new QLabel("No Cell");
    QL_noCell_3->setFont(QF_normal);
    QLabel *QL_type1_v_3 = new QLabel("Type 1");
    QL_type1_v_3->setFont(QF_normal);
    QLabel *QL_type2_v_3 = new QLabel("Type 2");
    QL_type2_v_3->setFont(QF_normal);
    QLabel *QL_type3_v_3 = new QLabel("Type 3");
    QL_type3_v_3->setFont(QF_normal);
    QLabel *QL_type1_h_3 = new QLabel("Type 1");
    QL_type1_h_3->setFont(QF_normal);
    QLabel *QL_type2_h_3 = new QLabel("Type 2");
    QL_type2_h_3->setFont(QF_normal);
    QLabel *QL_type3_h_3 = new QLabel("Type 3");
    QL_type3_h_3->setFont(QF_normal);
    
    // fourth matrix: properties of the cells
    QLabel *QL_K = new QLabel("K");
    QL_K->setFont(QF_normal);
    QLabel *QL_T_a = new QLabel("T_a");
    QL_T_a->setFont(QF_normal);
    QLabel *QL_T_b = new QLabel("T_b");
    QL_T_b->setFont(QF_normal);
    QLabel *QL_V0 = new QLabel("V0");
    QL_V0->setFont(QF_normal);
    QLabel *QL_type1 = new QLabel("Type 1");
    QL_type1->setFont(QF_normal);
    QLabel *QL_type2 = new QLabel("Type 2");
    QL_type2->setFont(QF_normal);
    QLabel *QL_type3 = new QLabel("Type 3");
    QL_type3->setFont(QF_normal);
    
    // create the content of the tables
    SpecialQDoubleSpinBox *QS_01_1 = new SpecialQDoubleSpinBox();
    QS_01_1->setMinimum(-10000);
    QS_01_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_02_1 = new SpecialQDoubleSpinBox();
    QS_02_1->setMinimum(-10000);
    QS_02_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_03_1 = new SpecialQDoubleSpinBox();
    QS_03_1->setMinimum(-10000);
    QS_03_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_11_1 = new SpecialQDoubleSpinBox();
    QS_11_1->setMinimum(-10000);
    QS_11_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_12_1 = new SpecialQDoubleSpinBox();
    QS_12_1->setMinimum(-10000);
    QS_12_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_13_1 = new SpecialQDoubleSpinBox();
    QS_13_1->setMinimum(-10000);
    QS_13_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_22_1 = new SpecialQDoubleSpinBox();
    QS_22_1->setMinimum(-10000);
    QS_22_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_23_1 = new SpecialQDoubleSpinBox();
    QS_23_1->setMinimum(-10000);
    QS_23_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_33_1 = new SpecialQDoubleSpinBox();
    QS_33_1->setMinimum(-10000);
    QS_33_1->setMaximum(10000);
    
    SpecialQDoubleSpinBox *QS_01_2 = new SpecialQDoubleSpinBox();
    QS_01_2->setMinimum(-10000);
    QS_01_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_02_2 = new SpecialQDoubleSpinBox();
    QS_02_2->setMinimum(-10000);
    QS_02_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_03_2 = new SpecialQDoubleSpinBox();
    QS_03_2->setMinimum(-10000);
    QS_03_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_11_2 = new SpecialQDoubleSpinBox();
    QS_11_2->setMinimum(-10000);
    QS_11_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_12_2 = new SpecialQDoubleSpinBox();
    QS_12_2->setMinimum(-10000);
    QS_12_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_13_2 = new SpecialQDoubleSpinBox();
    QS_13_2->setMinimum(-10000);
    QS_13_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_22_2 = new SpecialQDoubleSpinBox();
    QS_22_2->setMinimum(-10000);
    QS_22_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_23_2 = new SpecialQDoubleSpinBox();
    QS_23_2->setMinimum(-10000);
    QS_23_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_33_2 = new SpecialQDoubleSpinBox();
    QS_33_2->setMinimum(-10000);
    QS_33_2->setMaximum(10000);
    
    SpecialQDoubleSpinBox *QS_01_3 = new SpecialQDoubleSpinBox();
    QS_01_3->setMinimum(-10000);
    QS_01_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_02_3 = new SpecialQDoubleSpinBox();
    QS_02_3->setMinimum(-10000);
    QS_02_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_03_3 = new SpecialQDoubleSpinBox();
    QS_03_3->setMinimum(-10000);
    QS_03_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_11_3 = new SpecialQDoubleSpinBox();
    QS_11_3->setMinimum(-10000);
    QS_11_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_12_3 = new SpecialQDoubleSpinBox();
    QS_12_3->setMinimum(-10000);
    QS_12_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_13_3 = new SpecialQDoubleSpinBox();
    QS_13_3->setMinimum(-10000);
    QS_13_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_22_3 = new SpecialQDoubleSpinBox();
    QS_22_3->setMinimum(-10000);
    QS_22_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_23_3 = new SpecialQDoubleSpinBox();
    QS_23_3->setMinimum(-10000);
    QS_23_3->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_33_3 = new SpecialQDoubleSpinBox();
    QS_33_3->setMinimum(-10000);
    QS_33_3->setMaximum(10000);
    
    // Bulk Modulus K
    SpecialQDoubleSpinBox *QS_K_1 = new SpecialQDoubleSpinBox();
    QS_K_1->setMinimum(0);
    QS_K_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_K_2 = new SpecialQDoubleSpinBox();
    QS_K_2->setMinimum(0);
    QS_K_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_K_3 = new SpecialQDoubleSpinBox();
    QS_K_3->setMinimum(0);
    QS_K_3->setMaximum(10000);
    
    // Apical Surface Tension T_a
    SpecialQDoubleSpinBox *QS_T_a_1 = new SpecialQDoubleSpinBox();
    QS_T_a_1->setMinimum(-10000);
    QS_T_a_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_T_a_2 = new SpecialQDoubleSpinBox();
    QS_T_a_2->setMinimum(-10000);
    QS_T_a_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_T_a_3 = new SpecialQDoubleSpinBox();
    QS_T_a_3->setMinimum(-10000);
    QS_T_a_3->setMaximum(10000);

    // Basal Surface Tension T_b
    SpecialQDoubleSpinBox *QS_T_b_1 = new SpecialQDoubleSpinBox();
    QS_T_b_1->setMinimum(-10000);
    QS_T_b_1->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_T_b_2 = new SpecialQDoubleSpinBox();
    QS_T_b_2->setMinimum(-10000);
    QS_T_b_2->setMaximum(10000);
    SpecialQDoubleSpinBox *QS_T_b_3 = new SpecialQDoubleSpinBox();
    QS_T_b_3->setMinimum(-10000);
    QS_T_b_3->setMaximum(10000);
   
    // Preferred Volume V0
    SpecialQDoubleSpinBox *QS_V0_1 = new SpecialQDoubleSpinBox();
    QS_V0_1->setMinimum(0);
    QS_V0_1->setMaximum(1000000);
    SpecialQDoubleSpinBox *QS_V0_2 = new SpecialQDoubleSpinBox();
    QS_V0_2->setMinimum(0);
    QS_V0_2->setMaximum(1000000);
    SpecialQDoubleSpinBox *QS_V0_3 = new SpecialQDoubleSpinBox();
    QS_V0_3->setMinimum(0);
    QS_V0_3->setMaximum(1000000);
    
    QPushButton *QP_apply= new QPushButton("Apply changes");
    QGridLayout *layout = new QGridLayout;
	
    
    layout->setColumnStretch(6, 1);
    layout->setRowStretch(3, 1);
   
    // T_l matrix
    layout->addWidget(QL_T_l, 0, 2);
    layout->addWidget(QL_noCell_1, 2, 0);
    layout->addWidget(QL_type1_v_1, 3, 0);
    layout->addWidget(QL_type2_v_1, 4, 0);
    layout->addWidget(QL_type3_v_1, 5, 0);
    layout->addWidget(QL_type1_h_1, 1, 1);
    layout->addWidget(QL_type2_h_1, 1, 2);
    layout->addWidget(QL_type3_h_1, 1, 3);
    
    layout->addWidget(QS_01_1, 2, 1);
    layout->addWidget(QS_02_1, 2, 2);
    layout->addWidget(QS_03_1, 2, 3);
    layout->addWidget(QS_11_1, 3, 1);
    layout->addWidget(QS_12_1, 3, 2);
    layout->addWidget(QS_13_1, 3, 3);
    layout->addWidget(QS_22_1, 4, 2);
    layout->addWidget(QS_23_1, 4, 3);
    layout->addWidget(QS_33_1, 5, 3);
    
    // G_a matrix
    layout->addWidget(QL_G_a, 0, 9);
    layout->addWidget(QL_noCell_2, 2, 7);
    layout->addWidget(QL_type1_v_2, 3, 7);
    layout->addWidget(QL_type2_v_2, 4, 7);
    layout->addWidget(QL_type3_v_2, 5, 7);
    layout->addWidget(QL_type1_h_2, 1, 8);
    layout->addWidget(QL_type2_h_2, 1, 9);
    layout->addWidget(QL_type3_h_2, 1, 10);
    
    layout->addWidget(QS_01_2, 2, 8);
    layout->addWidget(QS_02_2, 2, 9);
    layout->addWidget(QS_03_2, 2, 10);
    layout->addWidget(QS_11_2, 3, 8);
    layout->addWidget(QS_12_2, 3, 9);
    layout->addWidget(QS_13_2, 3, 10);
    layout->addWidget(QS_22_2, 4, 9);
    layout->addWidget(QS_23_2, 4, 10);
    layout->addWidget(QS_33_2, 5, 10);
    
    // G_b matrix
    layout->addWidget(QL_G_b, 0, 16);
    layout->addWidget(QL_noCell_3, 2, 14);
    layout->addWidget(QL_type1_v_3, 3, 14);
    layout->addWidget(QL_type2_v_3, 4, 14);
    layout->addWidget(QL_type3_v_3, 5, 14);
    layout->addWidget(QL_type1_h_3, 1, 15);
    layout->addWidget(QL_type2_h_3, 1, 16);
    layout->addWidget(QL_type3_h_3, 1, 17);
    
    layout->addWidget(QS_01_3, 2, 15);
    layout->addWidget(QS_02_3, 2, 16);
    layout->addWidget(QS_03_3, 2, 17);
    layout->addWidget(QS_11_3, 3, 15);
    layout->addWidget(QS_12_3, 3, 16);
    layout->addWidget(QS_13_3, 3, 17);
    layout->addWidget(QS_22_3, 4, 16);
    layout->addWidget(QS_23_3, 4, 17);
    layout->addWidget(QS_33_3, 5, 17);
    
    // properties of the cells
    layout->addWidget(QL_K,     8, 0);
    layout->addWidget(QL_T_a,   9, 0);
    layout->addWidget(QL_T_b,  10, 0);
    layout->addWidget(QL_V0,   11, 0);
    layout->addWidget(QL_type1, 7, 1);
    layout->addWidget(QL_type2, 7, 2);
    layout->addWidget(QL_type3, 7, 3);

    layout->addWidget(QS_K_1,    8, 1);
    layout->addWidget(QS_K_2,    8, 2);
    layout->addWidget(QS_K_3,    8, 3);
    layout->addWidget(QS_T_a_1,  9, 1);
    layout->addWidget(QS_T_a_2,  9, 2);
    layout->addWidget(QS_T_a_3,  9, 3);
    layout->addWidget(QS_T_b_1, 10, 1);
    layout->addWidget(QS_T_b_2, 10, 2);
    layout->addWidget(QS_T_b_3, 10, 3);
    layout->addWidget(QS_V0_1,  11, 1);
    layout->addWidget(QS_V0_2,  11, 2);
    layout->addWidget(QS_V0_3,  11, 3);
    
    
    setLayout(layout);
    
    // connect the signals and slots
	//connect(QC_isPeriodic,SIGNAL(stateChanged(int)),this,SLOT(slot_isPeriodic(double)));
	
    raise();
    show();
    
    QApplication::processEvents();
}


void MechPropWindow::slot_show()
{
    if(isVisible()) {hide();}
    else {raise();show();}
}

//void MechPropWindow::signal_apply(double[4][3])
//{
//    
//    
//}



#endif // _USE_QT_
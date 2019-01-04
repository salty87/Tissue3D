/****************************************************************************
 **
 ** Copyright (C) 2008 Nokia Corporation and/or its subsidiary(-ies).
 ** Contact: Qt Software Information (qt-info@nokia.com)
 **
 ** This file is part of the example classes of the Qt Toolkit.
 **
 ** Commercial Usage
 ** Licensees holding valid Qt Commercial licenses may use this file in
 ** accordance with the Qt Commercial License Agreement provided with the
 ** Software or, alternatively, in accordance with the terms contained in
 ** a written agreement between you and Nokia.
 **
 **
 ** GNU General Public License Usage
 ** Alternatively, this file may be used under the terms of the GNU
 ** General Public License versions 2.0 or 3.0 as published by the Free
 ** Software Foundation and appearing in the file LICENSE.GPL included in
 ** the packaging of this file.  Please review the following information
 ** to ensure GNU General Public Licensing requirements will be met:
 ** http://www.fsf.org/licensing/licenses/info/GPLv2.html and
 ** http://www.gnu.org/copyleft/gpl.html.  In addition, as a special
 ** exception, Nokia gives you certain additional rights. These rights
 ** are described in the Nokia Qt GPL Exception version 1.3, included in
 ** the file GPL_EXCEPTION.txt in this package.
 **
 ** Qt for Windows(R) Licensees
 ** As a special exception, Nokia, as the sole copyright holder for Qt
 ** Designer, grants users of the Qt/Eclipse Integration plug-in the
 ** right for the Qt/Eclipse Integration to link to functionality
 ** provided by Qt Designer and its related libraries.
 **
 ** If you are unsure which license is appropriate for your use, please
 ** contact the sales department at qt-sales@nokia.com.
 **
 ****************************************************************************/

#define _STABILITY_CHECK_FLAT_TISSUE_

#ifdef _USE_QT_
  #include <QApplication>
  #include <QtGui>
  #include <QFont>
  #include <QGridLayout>
  #include <QPushButton>
  #include <QDialog>
  #include <QtCore/QCoreApplication>
  #include "TissueWindow.h"
  #include "MainWindow.h"
  #include "ProfileWindow.h"
  #include "HexagonsWindow.h"
#endif

#include "Integrator.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <QDir>
#include <QString>


int main(int argc, char *argv[])
{
    
//    Point P = Point(214,113,523);
//    
//    double ex = 200;
//    double ez = 1000;
//    
//    P.projectOnEllipse(ex, ez);
//    
//    return 1;
//
    //if(0){
    if(argc<4){
        std::cout << argc << std::endl;
        
        QApplication app(argc, argv);
        MainWindow widget;
        
        widget.activateWindow();
        
        srand(time(NULL));
        
        return app.exec();
    } else {
       // start Qt Engine in order to use QtScript
       QCoreApplication app(argc, argv);
        
       // now run simulations
       Integrator *integrator = new Integrator();
     
        //script name/Users/silvanus/MPI/EpithelialMechanics/SimulationResults/playground/run_1.3dvm_scr

        // iterate through all input files and run them
        for(int scripts_count = 1; scripts_count < argc-2; scripts_count++)
        {
           QString scriptFile(argv[scripts_count]);
//           QString scriptFile = "/Users/silvanus/MPI/EpithelialMechanics/SimulationResults/stabilityCheck/";
            
            std::cout << scriptFile.toStdString() << std::endl;
           
           // run the scriptz
           integrator->runScript(scriptFile);
        }
        
       delete integrator;
    }
}

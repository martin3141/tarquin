#ifndef TARQUIN_MAIN_WINDOW_INCLUDED 
#define TARQUIN_MAIN_WINDOW_INCLUDED 

#include "ui_mainwindow.h"
#include <qwt_plot.h>
#include <qwt_plot_picker.h>
#include <qwt_compat.h>
#include <QtWidgets>

#include "MRI.hpp"
#include "CFID.hpp"
#include <iostream>

class TarquinPlotWidget;
class Session;

class MyGraphicsScene : public QGraphicsScene
{
    Q_OBJECT

    public:
        MyGraphicsScene ( QWidget * parent = 0 );
        void SetVoxelFormat(const tarquin::coord& voxel, const QPen& pen, int Z);
};


class MyGraphicsView : public QGraphicsView
{
    Q_OBJECT

    public:
        MyGraphicsView ( QWidget * parent = 0 );
        MyGraphicsView ( QGraphicsScene * scene, QWidget * parent = 0 );

        //void SetVoxelFormat(const tarquin::coord& voxel, const QPen& pen, int Z);

    protected:
        virtual void wheelEvent ( QWheelEvent * event )
        {
            this->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

            // Scale the view / do the zoom
            double scaleFactor = 1.15;
            if(event->delta() > 0) {
                // Zoom in
                this->scale(scaleFactor, scaleFactor);
            } else {
                // Zooming out
                this->scale(1.0 / scaleFactor, 1.0 / scaleFactor);
            }
        }
        //void mousePressEvent ( QMouseEvent * event );

    public:
        QPen selected_pen;
        QPen in_fit_list_pen;
        QPen std_vox_pen;
        QBrush selected_brush;
        QBrush in_fit_list_brush;
        QBrush std_vox_brush;
};

typedef double coord_t;         // coordinate type
typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2
 
struct Point {
        coord_t x, y;
 
        bool operator <(const Point &p) const {
                return x < p.x || (x == p.x && y < p.y);
        }
};




class MainWindow : public QMainWindow
{
	Q_OBJECT

public:

	MainWindow(QWidget *parent = 0, Qt::WFlags flags = 0);

	TarquinPlotWidget* GetPlot() { return m_plot; }

private:

	Session* MakeNewSession();
	
    Session* KillSession();

	void SetWindowTitle(QString filename);

private slots:
	
	void OnFileQF();

	void OnFileOpenFID();

	void OnFileAbout();
	
	void OnFilePrint();
	
    void OnFilePrintLocalisation();
	
	void OnFileExportPDF();
	
    void OnFileExportJPEG();
    
    void OnFileExportPNG();
    
    void OnFileExportLocalisationPDF();
    
    void OnFileExportLocalisationPNG();
    
    void OnFileExportLocalisationJPEG();
    
    void OnFileExportWindowPNG();

	void OnFileExit();
	
	void OnAnalysisRunTARQUIN();
	
    void OnAnalysisPreprocess();

	void OnViewReal();
	
	void OnViewImag();
	
    void OnViewAbs();
	
	void OnViewTimeDomainSec();
	
    void OnViewTimeDomainPts();
	
    void OnViewFreqDomainPPM();
	
	void OnViewFreqDomainHz();
	
    void OnViewFreqDomainPts();

	void OnNextFit();
	
	void OnPreviousFit();

    void OnSubtractBaseline();
    
    void OnShowRawData();
    
    void OnSetLB();
    
    void OnSetZF();

    void OnSetBS();
    
    void OnAAToggle();

    void OnWUSToggle();

    void OnSetLeftppm();
    
    void OnSetRightppm();

	void OnViewAmplitudes();

	void OnExportPDF();
	
    void OnExportPDFStack();

	void OnExportTXT();

	void OnExportCSV();
	
    void OnExportCSV_SV();

	void OnExportCSVFit();
	
    void OnExportCSVFitMag();
	
    void OnExportCSVFit_SV();
	
    void OnExportCSVSpectra();
    
    void OnExportCSVSpectraMag();

    void OnExportDPT(bool processed);
    
    void OnExportDPTRaw();
    
    void OnExportDPTProc();
    
    void OnExportDPTWater();
	
    void OnAutoX();
    
    void OnAutoY();
    
    void OnPhi0();
    
    void OnPhi1();

    void OnVoxelChange(bool plot_update = true);

    void OnFitChange();
    
    void OnSelNone();
    
    void OnSelAll();

    void OnCloseMRI();

    void OnWWLL();
    
    void OnHideGrid();
    
    void OnHideLines();
    
    void OnHideMRI();
    
    void OnSetTrans();
    
    void OnLoadMRI();
    
    void UpdateGeom();

    void OnEditFin();

    void OnCenterVoxel();

    void mouseMoveEvent();

    //void wheelEvent(QWheelEvent* event);

    void mousePressEvent ( QMouseEvent * event );

    void OnCalPPM();
    
    void OnSetPivot();

    void PPMSelected(const QwtDoublePoint &point);

    void PivotSelected(const QwtDoublePoint &point);

    bool eventFilter(QObject *obj, QEvent *event);
    
    void OnListChange();

    void OnSliderChange();
    
    void OnSpinChange();
    
    QColor GetColMap(double val);

    private:

	bool run_preprocessor();

	bool run_tarquin();
	
	bool simulate_basis();

    std::vector<Point> convex_hull(std::vector<Point> P);

    coord2_t cross(const Point &O, const Point &A, const Point &B);

    Point calc_pt(double d, Eigen::Vector3d line_norm, Eigen::Vector3d line_point, Eigen::Matrix3d M_cut);
    
    // this version works better in some cases
    Point calc_pt(double d, Eigen::Vector3d line_norm, Eigen::Vector3d line_point, const Eigen::Matrix4d& M);
  

    typedef struct {
        double r,g,b;
    } COLOUR;
    

    COLOUR GetColour(double v,double vmin,double vmax)
    {
        COLOUR c = {1.0,1.0,1.0}; // white
        double dv;

        if (v < vmin)
            v = vmin;
        if (v > vmax)
            v = vmax;
        dv = vmax - vmin;

        if (v < (vmin + 0.25 * dv)) {
            c.r = 0;
            c.g = 4 * (v - vmin) / dv;
        } else if (v < (vmin + 0.5 * dv)) {
            c.r = 0;
            c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
        } else if (v < (vmin + 0.75 * dv)) {
            c.r = 4 * (v - vmin - 0.5 * dv) / dv;
            c.b = 0;
        } else {
            c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
            c.b = 0;
        }

        return(c);
    }

    double interpolate( double val, double y0, double x0, double y1, double x1 );

    double base( double val );

    double red( double gray );
    
    double green( double gray );
    
    double blue( double gray );

private:

	//! The user interface for this program.
	Ui::MainWindow m_ui;

	//! The plotting widget.
	TarquinPlotWidget* m_plot;
    
    //! MRS Goem view
    MyGraphicsView* m_view;
    
    //! MRS Goem scene
    MyGraphicsScene* m_scene;
    
	//! The current session.
	Session* m_session;

    //! Temporary mouse drag starting position
    int m_start_pos_y;
    int m_start_pos_x;

    //! Phi0 or Phi1
    bool m_phi0;
    bool m_phi1;
    
    //! Change WW/WL?
    bool m_ww_wl;

};



#endif // TARQUIN_MAIN_WINDOW_INCLUDED 

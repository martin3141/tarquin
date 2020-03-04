#ifndef TARQUIN_SESSION_INCLUDED
#define TARQUIN_SESSION_INCLUDED

// Q_MOC_RUN - workaround for a bug in Qt4 moc
// https://bugreports.qt.io/browse/QTBUG-22829
#ifndef Q_MOC_RUN
#include "Workspace.hpp"
#endif
#include <QObject>
#include <qwt_double_rect.h>
#include "MRI.hpp"

class MainWindow;

class QString;
class QColor;
class QwtPlotItem;

enum show_curve_type
{
	DISPLAY_CURVE,
	HIDE_CURVE
};


class Session : public QObject
{
	Q_OBJECT

	public:

		Session(MainWindow* parent);

		enum plot_domain_e
		{
			TIME_DOMAIN,
			FREQUENCY_DOMAIN
		};

		enum complex_mode_e
		{
			REAL,
			IMAG,
			ABS
		};

		enum units_fd_e
		{
			PPM,
			HERTZ,
			FD_PTS
		};

		enum units_td_e
		{
			SEC,
			TD_PTS
		};

		/*!
		 * Options that control the plotting, the user can control plotting through
		 * the plot widget, but these flags are more to do with what can be shown,
		 * i.e. what results are available for the user to toggle.
		 */
		struct ShowFlags
		{
			/*! time or frequency domain plot */
			plot_domain_e      domain;

			/*! plot real, imag or abs? */
			complex_mode_e     mode;

			/*! if frequency domain, what are we showing on the x axis */
			units_fd_e   units_fd;

			/*! if time domain, what are we showing on the x axis */
			units_td_e   units_td;

			/*! display the input signal */
			bool               show_input;

			/*! display the preprocessed signal */
			bool               show_preproc;

			/*! display the fitted signal */
			bool               show_model;
		    
            /*! subtract the baseline */
            bool               SB;

            /*! anti-aliasaing */
            bool               AA;
            
            /* plot WUS data instead? */
            bool               WUS;

			/*! indices of basis vectors to display */
			std::vector<int>   basis_indices;

            /* Phi1 pivot point in Hz */
            double             pivot;

            ShowFlags() : 
                domain(TIME_DOMAIN),
                mode(REAL),
                units_fd(HERTZ),
                units_td(SEC),
                show_input(false), 
                show_preproc(false),
                show_model(false),
                SB(false),
                AA(true),
                WUS(false),
                pivot(0.0)
            { }

        };

		void Update(bool all=true);

		void ShowResults();

		bool FitAvailable();

		void Save(const QString& filename);

		void Load(const QString& filename);
		
        void SetVoxel(const tarquin::coord& voxel);

        const tarquin::coord GetVoxel();

	public:

		tarquin::Workspace& GetWorkspace() { return m_workspace; }

        cvm::rvector map;

		//! The flags controlling what is shown.
		ShowFlags m_show_flags;
		
		//! fit number to be displayed
		int	m_fit_number;

        //! current display voxel
        tarquin::coord m_voxel;

        //! is the current display voxel in the fit list?
        bool m_in_fit_list;

        //! state of the processing
        bool data_loaded;
        bool data_preprocessed;
        bool data_fitted;
        
        bool mri_loaded;
        int grid_trans;

        MRI image;

	private:

		void add_timedomain_signal(
				const QString&           label, 
				const QColor&            colour,
				const float              dt,    
				const int                ns,    
				const int                ne,    
				const cvm::cvector&      y,
				show_curve_type          show_curve=DISPLAY_CURVE
				);

		void add_frequencydomain_signal(
				const QString&           label, 
				const QColor&            colour,
				cvm::rvector&            hz,
				cvm::rvector&            ppm,
				const cvm::cvector&      y,
				show_curve_type          show_curve=DISPLAY_CURVE
				);

	private slots:
		

	private:

		//! The main window.
		MainWindow* m_parent;

		//! The workspace.
		tarquin::Workspace m_workspace;

};

#endif // TARQUIN_SESSION_INCLUDED

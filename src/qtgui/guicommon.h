#ifndef TARQUIN_GUI_COMMON_INCLUDED
#define TARQUIN_GUI_COMMON_INCLUDED

class QString;
class QWidget;

namespace tarquin
{
class CBoswell;
} // namespace tarquin

void ErrorDialog(QWidget* parent, const QString& title, const QString& message);

void InfoDialog(QWidget* parent, const QString& title, const QString& message);

tarquin::CBoswell& GetLog();

#endif // TARQUIN_GUI_COMMON_INCLUDED

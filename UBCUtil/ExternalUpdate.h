#ifndef EXTERNAL_UPDATE_H
#define EXTERNAL_UPDATE_H

class ExternalUpdate {
public:
	virtual void update() = 0;
	virtual bool hasFutureUpdates() = 0;

	virtual void onStartup() = 0;
};

#endif
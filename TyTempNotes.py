log = open('fahtemp.log','r')
df = pd.read_csv(log)
df.columns = ["date","temp","target","status"]
df['date'] = pd.to_datetime(df['date'],unit='s')

run = pd.DataFrame(df[df['status'] == ' FAH is running']['temp'])
notr = pd.DataFrame(df[df['status'] == ' FAH is not running']['temp'])

run=run.rolling(21,min_periods=1).mean()
notr=notr.rolling(21,min_periods=1).mean()
temp7=pd.DataFrame(temp.rolling(21,min_periods=1).mean())

df['temp']=df['temp'].rolling(21,min_periods=1).mean()

run.insert(loc=0,column='date',value=df.take(run.index)['date'])
notr.insert(loc=0,column='date',value=df.take(notr.index)['date'])
temp7.insert(loc=0,column='date',value=df.take(temp7.index)['date'])

fig = plt.figure()
axes = fig.add_subplot(111)
axes.set_ylim(60,75)
plt.plot(run['temp'],color='r')
plt.plot(notr['temp'],color='b')

plt.scatter(run['date'], run['temp'],color='r')
plt.scatter(notr['date'],notr['temp'],color='b')

#plt.plot(temp7)
plt.show()



#plt.plot(df.iloc[:,1])
#plt.plot(temp7)

#ax = plt.subplots()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(60,75)
colors = {' FAH is running':'r', ' FAH is not running':'b'}
for n, g in df.groupby('status'):
     #g.plot('date','temp', ax=ax, color=colors[n])
     g.plot.scatter('date','temp', marker='.',ax=ax, color=colors[n])
#plt.scatter(temp7['date'],temp7['temp'],c='black')
plt.tight_layout()     
plt.show() 


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(60,75)

for n in range(len(df['status'])):
    if df['status'].loc[n]==" FAH is running":
        plt.plot(df['date'].loc[n],df['temp'].loc[n],c='r')
    elif df['status'].loc[n]==" FAH is not running":
        plt.plot(df['date'].loc[n],df['temp'].loc[n],c='b')

plt.tight_layout()     
plt.show() 


### WORKING changing color shit
ts = [0]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim(60,75)
i=0
for n in range(len(df['status'])-1):
    if df['status'].loc[n+1]==" FAH is running" and df['status'].loc[n]==" FAH is not running":
        ts.append(n)
        x = [df['temp'][ts[i]],df['temp'][n]]
        y = [df['date'][ts[i]],df['date'][n]]
        i+=1
        plt.plot(y,x,c='r')
        
    elif df['status'].loc[n+1]==" FAH is not running" and df['status'].loc[n]==" FAH is running":
        ts.append(n)
        x = [df['temp'][ts[i]],df['temp'][n]]
        y = [df['date'][ts[i]],df['date'][n]]
        i+=1
        plt.plot(y,x,c='b')

plt.plot(df['date'],df['target'],color='black')
plt.tight_layout()     
plt.show() 
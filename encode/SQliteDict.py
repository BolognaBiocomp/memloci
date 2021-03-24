from sqlite3 import dbapi2 as sqlite
import UserDict,  os
import pickle

class SQliteDict(UserDict.DictMixin):
    ''' dbdict, a dictionnary-like object for large datasets (several Tera-bytes) '''
    def __init__(self,dictName):
        self.db_filename = "%s.sqlite" % dictName
        if not os.path.isfile(self.db_filename):
            self.con = sqlite.connect(self.db_filename)
            self.con.execute("create table data (key PRIMARY KEY,value)")
        else:
            self.con = sqlite.connect(self.db_filename)
    def __getitem__(self, key):
        row = self.con.execute("select value from data where key=?",(key,)).fetchone()
        if not row: raise KeyError(key)
        return pickle.loads(str(row[0]))
    def __setitem__(self, key, item):
        if self.con.execute("select key from data where key=?",(key,)).fetchone():
            self.con.execute("update data set value=? where key=?",(pickle.dumps(item),key))
        else:
            self.con.execute("insert into data (key,value) values (?,?)",(key, pickle.dumps(item)))
        write_errors=0
        while write_errors<=5:
            try:
                self.con.commit()
                break
            except:
                os.wait(1)
                write_errors+=1
        if write_errors>5:
            self.con.commit()
    def __delitem__(self, key):
        if self.con.execute("select key from data where key=?",(key,)).fetchone():
            self.con.execute("delete from data where key=?",(key,))
            write_errors=0
            while write_errors<=5:
                try:
                    self.con.commit()
                    break
                except:
                    os.wait(1)
                    write_errors+=1
            if write_errors>5:
                self.con.commit()
        else:   raise KeyError(key)
    def keys(self):
        return [row[0] for row in self.con.execute("select key from data").fetchall()]

"""Measurements database for Quantum Explorer.

Usage in experiment scripts:
    from db_utils import record, query, query_table

    # After computing a result:
    record(sprint=65, model='hybrid', q=5, n=8, quantity='c',
           value=1.10, error=0.10, method='entropy_profile')

    # Look up prior results:
    rows = query(quantity='c', q=5)
    for r in rows:
        print(r)  # (id, sprint, model, q, n, quantity, value, error, method, notes)

    # Pretty-print a comparison table:
    query_table(quantity='c', model='hybrid')
"""
import sqlite3
import os

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results.db')

def _get_conn():
    conn = sqlite3.connect(DB_PATH)
    conn.execute('''CREATE TABLE IF NOT EXISTS measurements (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        sprint INTEGER,
        model TEXT,
        q INTEGER,
        n INTEGER,
        quantity TEXT,
        value REAL,
        error REAL,
        method TEXT,
        notes TEXT
    )''')
    return conn

def record(sprint, model, q, n, quantity, value, error=None, method=None, notes=None):
    """Record a measurement. Upserts: same (sprint, model, q, n, quantity) replaces old value."""
    conn = _get_conn()
    # Delete any existing row with same key
    conn.execute(
        'DELETE FROM measurements WHERE sprint=? AND model=? AND q=? AND n=? AND quantity=?',
        (sprint, model, q, n, quantity))
    conn.execute(
        'INSERT INTO measurements (sprint, model, q, n, quantity, value, error, method, notes) '
        'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
        (sprint, model, q, n, quantity, value, error, method, notes))
    conn.commit()
    conn.close()

def query(quantity=None, q=None, model=None, n=None, sprint=None):
    """Query measurements. Returns list of tuples (id, sprint, model, q, n, quantity, value, error, method, notes).
    All parameters are optional filters."""
    conn = _get_conn()
    sql = 'SELECT * FROM measurements WHERE 1=1'
    params = []
    if quantity is not None: sql += ' AND quantity=?'; params.append(quantity)
    if q is not None: sql += ' AND q=?'; params.append(q)
    if model is not None: sql += ' AND model=?'; params.append(model)
    if n is not None: sql += ' AND n=?'; params.append(n)
    if sprint is not None: sql += ' AND sprint=?'; params.append(sprint)
    sql += ' ORDER BY q, sprint DESC'
    rows = conn.execute(sql, params).fetchall()
    conn.close()
    return rows

def query_table(quantity=None, model=None, q=None):
    """Print a formatted table of measurements."""
    rows = query(quantity=quantity, model=model, q=q)
    if not rows:
        print("No matching measurements.")
        return
    print(f"{'sprint':>6} {'model':>8} {'q':>3} {'n':>4} {'quantity':>12} {'value':>12} {'error':>8} {'method'}")
    print("-" * 75)
    for r in rows:
        _, sprint, mdl, q_val, n_val, qty, val, err, meth, notes = r
        err_str = f"{err:.4f}" if err else "—"
        meth_str = meth or "—"
        print(f"{sprint:>6} {mdl:>8} {q_val:>3} {n_val:>4} {qty:>12} {val:>12.5f} {err_str:>8} {meth_str}")

def best(quantity, q, model='hybrid'):
    """Get the best (most recent, largest n) measurement of a quantity."""
    conn = _get_conn()
    row = conn.execute(
        'SELECT value, error, n, sprint, method FROM measurements '
        'WHERE quantity=? AND q=? AND model=? ORDER BY n DESC, sprint DESC LIMIT 1',
        (quantity, q, model)).fetchone()
    conn.close()
    if row:
        return {'value': row[0], 'error': row[1], 'n': row[2], 'sprint': row[3], 'method': row[4]}
    return None

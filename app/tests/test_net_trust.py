"""Tests for the OS-trust-store helper (app/net_trust.py)."""
import ssl
import sys

import pytest

sys.path.insert(0, __file__.rsplit("/tests/", 1)[0])
from net_trust import enable_os_trust_store  # noqa: E402


def test_enable_returns_tuple_and_never_raises():
    ok, msg = enable_os_trust_store()
    assert isinstance(ok, bool)
    assert isinstance(msg, str) and msg


def test_enabled_when_truststore_present():
    pytest.importorskip("truststore")
    ok, msg = enable_os_trust_store()
    assert ok is True
    assert "OS trust store" in msg


def test_verification_is_never_disabled():
    """The whole point: we change WHERE roots come from, never whether we verify."""
    enable_os_trust_store()
    ctx = ssl.create_default_context()
    assert ctx.verify_mode == ssl.CERT_REQUIRED
    assert ctx.check_hostname is True


def test_missing_truststore_degrades_safely(monkeypatch):
    """If truststore isn't installed we must no-op, not crash and not weaken TLS."""
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "truststore":
            raise ImportError("simulated missing truststore")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    ok, msg = enable_os_trust_store()
    assert ok is False and "certifi" in msg
    ctx = ssl.create_default_context()
    assert ctx.verify_mode == ssl.CERT_REQUIRED  # still verifying


def test_inject_failure_degrades_safely(monkeypatch):
    """A truststore that imports but fails to inject must not break startup."""
    truststore = pytest.importorskip("truststore")
    monkeypatch.setattr(truststore, "inject_into_ssl",
                        lambda: (_ for _ in ()).throw(RuntimeError("boom")))
    ok, msg = enable_os_trust_store()
    assert ok is False and "could not be enabled" in msg

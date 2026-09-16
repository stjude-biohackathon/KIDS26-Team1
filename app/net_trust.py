"""Use the operating-system trust store for TLS instead of certifi's bundle.

Why this exists
---------------
Corporate TLS-inspection proxies (Cloudflare Zero Trust / WARP, Zscaler,
Netskope, many enterprise VPNs) re-sign HTTPS traffic with a private root CA.
That CA is normally installed in the OS trust store (macOS Keychain, Windows
cert store), so `curl` and browsers work fine -- but Python's `ssl` module
defaults to the `certifi` bundle, which does NOT contain it. The symptom is a
host-specific failure like:

    URLError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed:
    self-signed certificate in certificate chain

...affecting both urllib and requests, while other hosts work normally.

`truststore` redirects Python's TLS verification to the OS trust store, so the
proxy's CA is honoured exactly as the browser honours it.

IMPORTANT: this does NOT disable or weaken certificate verification. Verification
stays fully on; only the source of trusted roots changes. We deliberately never
fall back to an unverified context -- if truststore is unavailable we leave
Python's default behaviour untouched and let the caller surface a real error.
"""
from __future__ import annotations


def enable_os_trust_store() -> tuple[bool, str]:
    """Best-effort: route TLS verification through the OS trust store.

    Returns (enabled, message). Never raises, and never weakens TLS -- a failure
    just leaves the stock certifi-based behaviour in place.
    """
    try:
        import truststore
    except ImportError:
        return False, ("truststore not installed; using certifi. If you are behind a "
                       "TLS-inspecting proxy, HTTPS to some hosts may fail.")
    try:
        truststore.inject_into_ssl()
        return True, "TLS verification using the OS trust store (truststore)."
    except Exception as e:  # noqa: BLE001 - must never break app startup
        return False, f"truststore present but could not be enabled ({e}); using certifi."

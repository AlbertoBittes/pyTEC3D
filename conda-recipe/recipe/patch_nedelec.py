#!/usr/bin/env python3
"""
Apply NEDELEC_ONE support for PRISM15/PRISM18 and PYRAMID13/PYRAMID14
to an upstream libMesh checkout by editing source files directly.

- No PRISM6 / PYRAMID5 support is added.
- Uses elem->positive_edge_orientation(i)
- Prism shape_deriv uses switch(j) style.

Usage:
  python3 apply_nedelec_prism_pyramid.py /path/to/libmesh
"""

from __future__ import annotations
import sys
from pathlib import Path


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def backup_file(path: Path) -> None:
    bak = path.with_suffix(path.suffix + ".bak")
    if not bak.exists():
        bak.write_bytes(path.read_bytes())


def ensure_contains_once(text: str, needle: str, where: str) -> None:
    if needle not in text:
        die(f"Expected to find marker '{needle}' in {where}")


def insert_after_marker(text: str, marker: str, insertion: str, where: str) -> str:
    idx = text.find(marker)
    if idx < 0:
        die(f"Could not find marker in {where}: {marker}")
    idx_end = idx + len(marker)
    return text[:idx_end] + insertion + text[idx_end:]


def insert_before_marker(text: str, marker: str, insertion: str, where: str) -> str:
    idx = text.find(marker)
    if idx < 0:
        die(f"Could not find marker in {where}: {marker}")
    return text[:idx] + insertion + text[idx:]


def already_done(text: str) -> bool:
    return ("case PRISM15" in text) or ("PRISM15:" in text)


def patch_fe_interface(root: Path) -> None:
    path = root / "src/fe/fe_interface.C"
    if not path.exists():
        die(f"File not found: {path}")
    backup_file(path)
    text = path.read_text(encoding="utf-8", errors="strict")

    # Idempotency check
    if "case PRISM15:" in text and "case PYRAMID13:" in text:
        print(f"[skip] {path} already contains PRISM/PYRAMID NEDELEC_ONE cases")
        return

    # We patch the NEDELEC_ONE switch(el_t) list. We insert after "case TET10:" if present,
    # otherwise after "case QUAD9:" (older trees).
    if "case TET10:" in text:
        marker = "case TET10:\n"
    else:
        marker = "case QUAD9:\n"

    ensure_contains_once(text, marker, str(path))

    insertion = (
        "        case PRISM15:\n"
        "        case PRISM18:\n"
        "        case PYRAMID13:\n"
        "        case PYRAMID14:\n"
    )
    text2 = insert_after_marker(text, marker, insertion, str(path))
    path.write_text(text2, encoding="utf-8")
    print(f"[ok] patched {path}")


def patch_fe_nedelec_one(root: Path) -> None:
    path = root / "src/fe/fe_nedelec_one.C"
    if not path.exists():
        die(f"File not found: {path}")
    backup_file(path)
    text = path.read_text(encoding="utf-8", errors="strict")

    if "case PRISM15:" in text and "case PYRAMID13:" in text:
        print(f"[skip] {path} already contains PRISM/PYRAMID cases")
        return

    # --- Patch nedelec_one_n_dofs switch(t) ---
    # Insert before "case INVALID_ELEM:" inside nedelec_one_n_dofs
    ensure_contains_once(text, "unsigned int nedelec_one_n_dofs", str(path))
    ensure_contains_once(text, "case INVALID_ELEM:", str(path))

    dofs_insertion = (
        "\n"
        "        case PRISM15:\n"
        "        case PRISM18:\n"
        "          return 9;\n"
        "\n"
        "        case PYRAMID13:\n"
        "        case PYRAMID14:\n"
        "          return 8;\n"
    )
    # Only insert if not already there
    if "return 9;" not in text or "PRISM15" not in text:
        text = insert_before_marker(text, "        case INVALID_ELEM:", dofs_insertion, str(path))

    # --- Patch nedelec_one_n_dofs_at_node switch(t) ---
    ensure_contains_once(text, "unsigned int nedelec_one_n_dofs_at_node", str(path))

    # Insert before "case INVALID_ELEM:" inside the at_node switch(t).
    # We'll use a more specific anchor: the SECOND occurrence of "case INVALID_ELEM:" after the at_node function starts.
    start = text.find("unsigned int nedelec_one_n_dofs_at_node")
    if start < 0:
        die(f"Could not locate nedelec_one_n_dofs_at_node in {path}")
    inv1 = text.find("case INVALID_ELEM:", start)
    if inv1 < 0:
        die(f"Could not locate 'case INVALID_ELEM' in {path} within n_dofs_at_node")

    atnode_insertion = (
        "\n"
        "        case PRISM15:\n"
        "        {\n"
        "          switch (n)\n"
        "          {\n"
        "            case 0: case 1: case 2: case 3: case 4: case 5:\n"
        "              return 0;\n"
        "            case 6: case 7: case 8: case 9: case 10: case 11: case 12: case 13: case 14:\n"
        "              return 1;\n"
        "            default:\n"
        "              libmesh_error_msg(\"ERROR: Invalid node ID \" << n);\n"
        "          }\n"
        "        }\n"
        "\n"
        "        case PRISM18:\n"
        "        {\n"
        "          switch (n)\n"
        "          {\n"
        "            case 0: case 1: case 2: case 3: case 4: case 5:\n"
        "              return 0;\n"
        "            case 6: case 7: case 8: case 9: case 10: case 11: case 12: case 13: case 14:\n"
        "              return 1;\n"
        "            case 15: case 16: case 17:\n"
        "              return 0;\n"
        "            default:\n"
        "              libmesh_error_msg(\"ERROR: Invalid node ID \" << n);\n"
        "          }\n"
        "        }\n"
        "\n"
        "        case PYRAMID13:\n"
        "        {\n"
        "          switch (n)\n"
        "          {\n"
        "            case 0: case 1: case 2: case 3: case 4:\n"
        "              return 0;\n"
        "            case 5: case 6: case 7: case 8: case 9: case 10: case 11: case 12:\n"
        "              return 1;\n"
        "            default:\n"
        "              libmesh_error_msg(\"ERROR: Invalid node ID \" << n);\n"
        "          }\n"
        "        }\n"
        "\n"
        "        case PYRAMID14:\n"
        "        {\n"
        "          switch (n)\n"
        "          {\n"
        "            case 0: case 1: case 2: case 3: case 4:\n"
        "              return 0;\n"
        "            case 5: case 6: case 7: case 8: case 9: case 10: case 11: case 12:\n"
        "              return 1;\n"
        "            case 13:\n"
        "              return 0;\n"
        "            default:\n"
        "              libmesh_error_msg(\"ERROR: Invalid node ID \" << n);\n"
        "          }\n"
        "        }\n"
    )

    if "case PRISM15:" not in text:
        text = text[:inv1] + atnode_insertion + text[inv1:]

    path.write_text(text, encoding="utf-8")
    print(f"[ok] patched {path}")


def patch_fe_nedelec_one_shape_3D(root: Path) -> None:
    path = root / "src/fe/fe_nedelec_one_shape_3D.C"
    if not path.exists():
        die(f"File not found: {path}")
    backup_file(path)
    text = path.read_text(encoding="utf-8", errors="strict")

    if "case PRISM15:" in text and "case PYRAMID13:" in text and "shape_deriv" in text:
        if "case PRISM15:" in text and "case PRISM18:" in text and "case PYRAMID14:" in text:
            # Still could be partial; we proceed with targeted checks below
            pass

    # ---- Insert into shape() switch(elem->type()) before default ----
    ensure_contains_once(text, "RealGradient FE<3,NEDELEC_ONE>::shape", str(path))
    ensure_contains_once(text, "default:", str(path))

    # We'll add cases only if not already present
    if "case PRISM15:" not in text or "case PYRAMID13:" not in text:
        # Find the shape() block by locating its switch(elem->type())
        shape_start = text.find("RealGradient FE<3,NEDELEC_ONE>::shape")
        shape_switch = text.find("switch (elem->type())", shape_start)
        if shape_switch < 0:
            die(f"Could not find switch(elem->type()) in shape() in {path}")

        # Find the "default:" inside that switch (the first default after shape_switch)
        default_pos = text.find("\n        default:", shape_switch)
        if default_pos < 0:
            die(f"Could not find default: inside shape() switch in {path}")

        shape_cases = (
            "\n"
            "        case PRISM15:\n"
            "        case PRISM18:\n"
            "        {\n"
            "          libmesh_assert_less (i, 9);\n"
            "\n"
            "          const Real xi   = p(0);\n"
            "          const Real eta  = p(1);\n"
            "          const Real zeta = p(2);\n"
            "\n"
            "          // Reference prism: triangle in (xi,eta), zeta in [-1,1]\n"
            "          libmesh_assert_greater_equal (xi,  -10.*TOLERANCE);\n"
            "          libmesh_assert_greater_equal (eta, -10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal    (xi+eta, 1.0+10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal    (std::fabs(zeta), 1.0+10.*TOLERANCE);\n"
            "\n"
            "          const Real t0 = 1.0 - xi - eta;\n"
            "          const Real t1 = xi;\n"
            "          const Real t2 = eta;\n"
            "\n"
            "          const Real f0 = 0.5*(1.0 - zeta);\n"
            "          const Real f1 = 0.5*(1.0 + zeta);\n"
            "\n"
            "          RealGradient v;\n"
            "          switch (i)\n"
            "          {\n"
            "            case 0: v = RealGradient(f0*(1.0-eta), f0*xi, 0.0); break; // (0,1)\n"
            "            case 1: v = RealGradient(f0*(-eta),    f0*xi, 0.0); break; // (1,2)\n"
            "            case 2: v = RealGradient(f0*(eta),     f0*(1.0-xi), 0.0); break; // (0,2)\n"
            "            case 3: v = RealGradient(0.0, 0.0, 0.5*t0); break; // (0,3)\n"
            "            case 4: v = RealGradient(0.0, 0.0, 0.5*t1); break; // (1,4)\n"
            "            case 5: v = RealGradient(0.0, 0.0, 0.5*t2); break; // (2,5)\n"
            "            case 6: v = RealGradient(f1*(1.0-eta), f1*xi, 0.0); break; // (3,4)\n"
            "            case 7: v = RealGradient(f1*(-eta),    f1*xi, 0.0); break; // (4,5)\n"
            "            case 8: v = RealGradient(f1*(eta),     f1*(1.0-xi), 0.0); break; // (3,5)\n"
            "            default:\n"
            "              libmesh_error_msg(\"Invalid i = \" << i);\n"
            "          }\n"
            "\n"
            "          if (!elem->positive_edge_orientation(i))\n"
            "            v = -v;\n"
            "          return v;\n"
            "        }\n"
            "\n"
            "        case PYRAMID13:\n"
            "        case PYRAMID14:\n"
            "        {\n"
            "          libmesh_assert_less (i, 8);\n"
            "\n"
            "          const Real xi   = p(0);\n"
            "          const Real eta  = p(1);\n"
            "          const Real zeta = p(2);\n"
            "\n"
            "          // Reference pyramid: base zeta=0, apex zeta=1\n"
            "          libmesh_assert_less_equal (std::fabs(xi),  1.0+10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal (std::fabs(eta), 1.0+10.*TOLERANCE);\n"
            "          libmesh_assert_greater_equal (zeta, -10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal    (zeta,  1.0+10.*TOLERANCE);\n"
            "\n"
            "          Real a = 1.0 - zeta;\n"
            "          if (std::fabs(a) < 1e-14)\n"
            "            a = (a >= 0.0) ? 1e-14 : -1e-14;\n"
            "\n"
            "          const Real inva  = 1.0 / a;\n"
            "          const Real inva2 = inva * inva;\n"
            "\n"
            "          const Real l0 = 0.25*(a - xi - eta + xi*eta*inva);\n"
            "          const Real l1 = 0.25*(a + xi - eta - xi*eta*inva);\n"
            "          const Real l2 = 0.25*(a + xi + eta + xi*eta*inva);\n"
            "          const Real l3 = 0.25*(a - xi + eta - xi*eta*inva);\n"
            "          const Real l4 = zeta;\n"
            "\n"
            "          const RealGradient g0(0.25*(-1.0 + eta*inva), 0.25*(-1.0 + xi*inva), 0.25*(-1.0 + xi*eta*inva2));\n"
            "          const RealGradient g1(0.25*( 1.0 - eta*inva), 0.25*(-1.0 - xi*inva), 0.25*(-1.0 - xi*eta*inva2));\n"
            "          const RealGradient g2(0.25*( 1.0 + eta*inva), 0.25*( 1.0 + xi*inva), 0.25*(-1.0 + xi*eta*inva2));\n"
            "          const RealGradient g3(0.25*(-1.0 - eta*inva), 0.25*( 1.0 - xi*inva), 0.25*(-1.0 - xi*eta*inva2));\n"
            "          const RealGradient g4(0.0, 0.0, 1.0);\n"
            "\n"
            "          RealGradient v;\n"
            "          switch (i)\n"
            "          {\n"
            "            case 0: v = l0*g1 - l1*g0; break; // (0,1)\n"
            "            case 1: v = l1*g2 - l2*g1; break; // (1,2)\n"
            "            case 2: v = l2*g3 - l3*g2; break; // (2,3)\n"
            "            case 3: v = l3*g0 - l0*g3; break; // (3,0)\n"
            "            case 4: v = l0*g4 - l4*g0; break; // (0,4)\n"
            "            case 5: v = l1*g4 - l4*g1; break; // (1,4)\n"
            "            case 6: v = l2*g4 - l4*g2; break; // (2,4)\n"
            "            case 7: v = l3*g4 - l4*g3; break; // (3,4)\n"
            "            default:\n"
            "              libmesh_error_msg(\"Invalid i = \" << i);\n"
            "          }\n"
            "\n"
            "          if (!elem->positive_edge_orientation(i))\n"
            "            v = -v;\n"
            "          return v;\n"
            "        }\n"
        )

        text = text[:default_pos] + shape_cases + text[default_pos:]

    # ---- Insert into shape_deriv() switch(elem->type()) before default ----
    ensure_contains_once(text, "RealGradient FE<3,NEDELEC_ONE>::shape_deriv", str(path))
    deriv_start = text.find("RealGradient FE<3,NEDELEC_ONE>::shape_deriv")
    deriv_switch = text.find("switch (elem->type())", deriv_start)
    if deriv_switch < 0:
        die(f"Could not find switch(elem->type()) in shape_deriv() in {path}")
    deriv_default_pos = text.find("\n        default:", deriv_switch)
    if deriv_default_pos < 0:
        die(f"Could not find default: inside shape_deriv() switch in {path}")

    if "case PRISM15:" not in text[deriv_start:deriv_default_pos]:
        deriv_cases = (
            "\n"
            "        case PRISM15:\n"
            "        case PRISM18:\n"
            "        {\n"
            "          libmesh_assert_less (i, 9);\n"
            "\n"
            "          const Real xi   = p(0);\n"
            "          const Real eta  = p(1);\n"
            "          const Real zeta = p(2);\n"
            "\n"
            "          libmesh_assert_greater_equal (xi,  -10.*TOLERANCE);\n"
            "          libmesh_assert_greater_equal (eta, -10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal    (xi+eta, 1.0+10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal    (std::fabs(zeta), 1.0+10.*TOLERANCE);\n"
            "\n"
            "          const Real f0  = 0.5*(1.0 - zeta);\n"
            "          const Real f1  = 0.5*(1.0 + zeta);\n"
            "          const Real df0 = -0.5;\n"
            "          const Real df1 =  0.5;\n"
            "\n"
            "          RealGradient dv;\n"
            "          switch (i)\n"
            "          {\n"
            "            case 0:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0,  f0,  0.0); break;\n"
            "                case 1: dv = RealGradient(-f0, 0.0, 0.0); break;\n"
            "                default: dv = RealGradient(df0*(1.0-eta), df0*xi, 0.0); break;\n"
            "              }\n"
            "              break;\n"
            "            case 1:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0,  f0,  0.0); break;\n"
            "                case 1: dv = RealGradient(-f0, 0.0, 0.0); break;\n"
            "                default: dv = RealGradient(df0*(-eta), df0*xi, 0.0); break;\n"
            "              }\n"
            "              break;\n"
            "            case 2:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0, -f0, 0.0); break;\n"
            "                case 1: dv = RealGradient( f0, 0.0, 0.0); break;\n"
            "                default: dv = RealGradient(df0*eta, df0*(1.0-xi), 0.0); break;\n"
            "              }\n"
            "              break;\n"
            "            case 3:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0, 0.0, -0.5); break;\n"
            "                case 1: dv = RealGradient(0.0, 0.0, -0.5); break;\n"
            "                default: dv = RealGradient(); break;\n"
            "              }\n"
            "              break;\n"
            "            case 4:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0, 0.0, 0.5); break;\n"
            "                default: dv = RealGradient(); break;\n"
            "              }\n"
            "              break;\n"
            "            case 5:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 1: dv = RealGradient(0.0, 0.0, 0.5); break;\n"
            "                default: dv = RealGradient(); break;\n"
            "              }\n"
            "              break;\n"
            "            case 6:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0,  f1,  0.0); break;\n"
            "                case 1: dv = RealGradient(-f1, 0.0, 0.0); break;\n"
            "                default: dv = RealGradient(df1*(1.0-eta), df1*xi, 0.0); break;\n"
            "              }\n"
            "              break;\n"
            "            case 7:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0,  f1,  0.0); break;\n"
            "                case 1: dv = RealGradient(-f1, 0.0, 0.0); break;\n"
            "                default: dv = RealGradient(df1*(-eta), df1*xi, 0.0); break;\n"
            "              }\n"
            "              break;\n"
            "            case 8:\n"
            "              switch (j)\n"
            "              {\n"
            "                case 0: dv = RealGradient(0.0, -f1, 0.0); break;\n"
            "                case 1: dv = RealGradient( f1, 0.0, 0.0); break;\n"
            "                default: dv = RealGradient(df1*eta, df1*(1.0-xi), 0.0); break;\n"
            "              }\n"
            "              break;\n"
            "            default:\n"
            "              libmesh_error_msg(\"Invalid i = \" << i);\n"
            "          }\n"
            "\n"
            "          if (!elem->positive_edge_orientation(i))\n"
            "            dv = -dv;\n"
            "          return dv;\n"
            "        }\n"
            "\n"
            "        case PYRAMID13:\n"
            "        case PYRAMID14:\n"
            "        {\n"
            "          libmesh_assert_less (i, 8);\n"
            "\n"
            "          const Real xi   = p(0);\n"
            "          const Real eta  = p(1);\n"
            "          const Real zeta = p(2);\n"
            "\n"
            "          libmesh_assert_less_equal (std::fabs(xi),  1.0+10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal (std::fabs(eta), 1.0+10.*TOLERANCE);\n"
            "          libmesh_assert_greater_equal (zeta, -10.*TOLERANCE);\n"
            "          libmesh_assert_less_equal    (zeta,  1.0+10.*TOLERANCE);\n"
            "\n"
            "          Real a = 1.0 - zeta;\n"
            "          if (std::fabs(a) < 1e-14)\n"
            "            a = (a >= 0.0) ? 1e-14 : -1e-14;\n"
            "\n"
            "          const Real inva  = 1.0 / a;\n"
            "          const Real inva2 = inva * inva;\n"
            "          const Real inva3 = inva2 * inva;\n"
            "\n"
            "          const Real l0 = 0.25*(a - xi - eta + xi*eta*inva);\n"
            "          const Real l1 = 0.25*(a + xi - eta - xi*eta*inva);\n"
            "          const Real l2 = 0.25*(a + xi + eta + xi*eta*inva);\n"
            "          const Real l3 = 0.25*(a - xi + eta - xi*eta*inva);\n"
            "          const Real l4 = zeta;\n"
            "\n"
            "          const RealGradient g0(0.25*(-1.0 + eta*inva), 0.25*(-1.0 + xi*inva), 0.25*(-1.0 + xi*eta*inva2));\n"
            "          const RealGradient g1(0.25*( 1.0 - eta*inva), 0.25*(-1.0 - xi*inva), 0.25*(-1.0 - xi*eta*inva2));\n"
            "          const RealGradient g2(0.25*( 1.0 + eta*inva), 0.25*( 1.0 + xi*inva), 0.25*(-1.0 + xi*eta*inva2));\n"
            "          const RealGradient g3(0.25*(-1.0 - eta*inva), 0.25*( 1.0 - xi*inva), 0.25*(-1.0 - xi*eta*inva2));\n"
            "          const RealGradient g4(0.0, 0.0, 1.0);\n"
            "\n"
            "          RealGradient dg0, dg1, dg2, dg3, dg4;\n"
            "          dg4 = RealGradient();\n"
            "          if (j == 0)\n"
            "          {\n"
            "            dg0 = RealGradient(0.0,           0.25*inva,       0.25*eta*inva2);\n"
            "            dg1 = RealGradient(0.0,          -0.25*inva,      -0.25*eta*inva2);\n"
            "            dg2 = RealGradient(0.0,           0.25*inva,       0.25*eta*inva2);\n"
            "            dg3 = RealGradient(0.0,          -0.25*inva,      -0.25*eta*inva2);\n"
            "          }\n"
            "          else if (j == 1)\n"
            "          {\n"
            "            dg0 = RealGradient(0.25*inva,     0.0,             0.25*xi*inva2);\n"
            "            dg1 = RealGradient(-0.25*inva,    0.0,            -0.25*xi*inva2);\n"
            "            dg2 = RealGradient(0.25*inva,     0.0,             0.25*xi*inva2);\n"
            "            dg3 = RealGradient(-0.25*inva,    0.0,            -0.25*xi*inva2);\n"
            "          }\n"
            "          else\n"
            "          {\n"
            "            dg0 = RealGradient(0.25*eta*inva2, 0.25*xi*inva2,   0.5*xi*eta*inva3);\n"
            "            dg1 = RealGradient(-0.25*eta*inva2,-0.25*xi*inva2, -0.5*xi*eta*inva3);\n"
            "            dg2 = RealGradient(0.25*eta*inva2, 0.25*xi*inva2,   0.5*xi*eta*inva3);\n"
            "            dg3 = RealGradient(-0.25*eta*inva2,-0.25*xi*inva2, -0.5*xi*eta*inva3);\n"
            "          }\n"
            "\n"
            "          const auto deriv_whitney = [&](const Real la, const Real lb,\n"
            "                                         const RealGradient & ga, const RealGradient & gb,\n"
            "                                         const RealGradient & dga, const RealGradient & dgb) -> RealGradient\n"
            "          {\n"
            "            const Real dla = ga(j);\n"
            "            const Real dlb = gb(j);\n"
            "            return dla*gb + la*dgb - dlb*ga - lb*dga;\n"
            "          };\n"
            "\n"
            "          RealGradient dv;\n"
            "          switch (i)\n"
            "          {\n"
            "            case 0: dv = deriv_whitney(l0,l1,g0,g1,dg0,dg1); break;\n"
            "            case 1: dv = deriv_whitney(l1,l2,g1,g2,dg1,dg2); break;\n"
            "            case 2: dv = deriv_whitney(l2,l3,g2,g3,dg2,dg3); break;\n"
            "            case 3: dv = -(deriv_whitney(l0,l3,g0,g3,dg0,dg3)); break; // N = l3*g0 - l0*g3\n"
            "            case 4: dv = deriv_whitney(l0,l4,g0,g4,dg0,dg4); break;\n"
            "            case 5: dv = deriv_whitney(l1,l4,g1,g4,dg1,dg4); break;\n"
            "            case 6: dv = deriv_whitney(l2,l4,g2,g4,dg2,dg4); break;\n"
            "            case 7: dv = deriv_whitney(l3,l4,g3,g4,dg3,dg4); break;\n"
            "            default:\n"
            "              libmesh_error_msg(\"Invalid i = \" << i);\n"
            "          }\n"
            "\n"
            "          libmesh_ignore(l4);\n"
            "\n"
            "          if (!elem->positive_edge_orientation(i))\n"
            "            dv = -dv;\n"
            "          return dv;\n"
            "        }\n"
        )

        text = text[:deriv_default_pos] + deriv_cases + text[deriv_default_pos:]

    path.write_text(text, encoding="utf-8")
    print(f"[ok] patched {path}")


def main() -> None:
    if len(sys.argv) != 2:
        die("Usage: python3 apply_nedelec_prism_pyramid.py /path/to/libmesh")

    root = Path(sys.argv[1]).resolve()
    if not root.exists():
        die(f"Path does not exist: {root}")
    if not (root / "src").exists():
        die(f"Does not look like a libMesh repo root (missing src/): {root}")

    patch_fe_interface(root)
    patch_fe_nedelec_one(root)
    patch_fe_nedelec_one_shape_3D(root)

    print("\nAll done.")
    print("Backups written as *.bak next to modified files.")
    print("Now rebuild libMesh.")


if __name__ == "__main__":
    main()

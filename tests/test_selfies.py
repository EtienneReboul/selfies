import faulthandler
import random

import pytest
from rdkit.Chem import MolFromSmiles

import selfies as sf

faulthandler.enable()


@pytest.fixture()
def max_selfies_len():
    return 1000


@pytest.fixture()
def large_alphabet():
    alphabet = sf.get_semantic_robust_alphabet()
    alphabet.update([
        "[#Br]", "[#Branch1]", "[#Branch2]", "[#Branch3]", "[#C@@H1]",
        "[#C@@]", "[#C@H1]", "[#C@]", "[#C]", "[#Cl]", "[#F]", "[#H]", "[#I]",
        "[#NH1]", "[#N]", "[#O]", "[#P]", "[#Ring1]", "[#Ring2]", "[#Ring3]",
        "[#S]", "[/Br]", "[/C@@H1]", "[/C@@]", "[/C@H1]", "[/C@]", "[/C]",
        "[/Cl]", "[/F]", "[/H]", "[/I]", "[/NH1]", "[/N]", "[/O]", "[/P]",
        "[/S]", "[=Br]", "[=Branch1]", "[=Branch2]", "[=Branch3]", "[=C@@H1]",
        "[=C@@]", "[=C@H1]", "[=C@]", "[=C]", "[=Cl]", "[=F]", "[=H]", "[=I]",
        "[=NH1]", "[=N]", "[=O]", "[=P]", "[=Ring1]", "[=Ring2]", "[=Ring3]",
        "[=S]", "[Br]", "[Branch1]", "[Branch2]", "[Branch3]", "[C@@H1]",
        "[C@@]", "[C@H1]", "[C@]", "[C]", "[Cl]", "[F]", "[H]", "[I]", "[NH1]",
        "[N]", "[O]", "[P]", "[Ring1]", "[Ring2]", "[Ring3]", "[S]", "[\\Br]",
        "[\\C@@H1]", "[\\C@@]", "[\\C@H1]", "[\\C@]", "[\\C]", "[\\Cl]",
        "[\\F]", "[\\H]", "[\\I]", "[\\NH1]", "[\\N]", "[\\O]", "[\\P]",
        "[\\S]", "[nop]"
    ])
    return list(alphabet)


def test_random_selfies_decoder(trials, max_selfies_len, large_alphabet):
    """Tests that SELFIES that are generated by randomly stringing together
    symbols from the SELFIES alphabet are decoded into valid SMILES.
    """

    alphabet = tuple(large_alphabet)

    for _ in range(trials):

        # create random SELFIES and decode
        rand_len = random.randint(1, max_selfies_len)
        rand_selfies = "".join(random_choices(alphabet, k=rand_len))
        smiles = sf.decoder(rand_selfies)

        # check if SMILES is valid
        try:
            is_valid = MolFromSmiles(smiles, sanitize=True) is not None
        except Exception:
            is_valid = False

        err_msg = "SMILES: {}\n\t SELFIES: {}".format(smiles, rand_selfies)
        assert is_valid, err_msg


def test_nop_symbol_decoder(max_selfies_len, large_alphabet):
    """Tests that the '[nop]' symbol is always skipped over.
    """

    alphabet = list(large_alphabet)
    alphabet.remove("[nop]")

    for _ in range(100):

        # create random SELFIES with and without [nop]
        rand_len = random.randint(1, max_selfies_len)
        rand_mol = random_choices(alphabet, k=rand_len)
        rand_mol.extend(["[nop]"] * (max_selfies_len - rand_len))
        random.shuffle(rand_mol)

        with_nops = "".join(rand_mol)
        without_nops = with_nops.replace("[nop]", "")

        assert sf.decoder(with_nops) == sf.decoder(without_nops)


def test_get_semantic_constraints():
    constraints = sf.get_semantic_constraints()
    assert constraints is not sf.get_semantic_constraints()  # not alias
    assert "?" in constraints


def test_change_constraints_cache_clear():
    alphabet = sf.get_semantic_robust_alphabet()
    assert alphabet == sf.get_semantic_robust_alphabet()
    assert sf.decoder("[C][#C]") == "C#C"

    new_constraints = sf.get_semantic_constraints()
    new_constraints["C"] = 1
    sf.set_semantic_constraints(new_constraints)

    new_alphabet = sf.get_semantic_robust_alphabet()
    assert new_alphabet != alphabet
    assert sf.decoder("[C][#C]") == "CC"

    sf.set_semantic_constraints()  # re-set alphabet


def test_invalid_or_unsupported_smiles_encoder():
    malformed_smiles = [
        "",
        "(",
        "C(Cl)(Cl)CC[13C",
        "C(CCCOC",
        "C=(CCOC",
        "CCCC)",
        "C1CCCCC",
        "C(F)(F)(F)(F)(F)F",  # violates bond constraints
        "C=C1=CCCCCC1",  # violates bond constraints
        "CC*CC",  # uses wildcard
        "C$C",  # uses $ bond
        "S[As@TB1](F)(Cl)(Br)N",  # unrecognized chirality,
        "SOMETHINGWRONGHERE",
        "1243124124",
    ]

    for smiles in malformed_smiles:
        with pytest.raises(sf.EncoderError):
            sf.encoder(smiles)


def test_malformed_selfies_decoder():
    with pytest.raises(sf.DecoderError):
        sf.decoder("[O][=C][O][C][C][C][C][O][N][Branch2_3")


def random_choices(population, k):  # random.choices was new in Python v3.6
    return [random.choice(population) for _ in range(k)]


def test_decoder_attribution():
    sm, am = sf.decoder(
        "[C][N][C][Branch1][C][P][C][C][Ring1][=Branch1]", attribute=True)
    # check that P lined up
    for ta in am:
        if ta[0] == 'P':
            for i, v in ta[1]:
                if v == '[P]':
                    return
    raise ValueError('Failed to find P in attribution map')


def test_encoder_attribution():
    smiles = "C1([O-])C=CC=C1Cl"
    indices = [0, 3, 3, 3, 5, 7, 8, 10, None, None, 12]
    s, am = sf.encoder(smiles, attribute=True)
    # check that Cl lined up
    for i, ta in enumerate(am):
        if ta[1]:
            assert indices[i] == ta[1][0][0], f'found {ta[1]}; should be {indices[i]}'
        if ta[0] == '[Cl]':
            for i, v in ta[1]:
                if v == 'Cl':
                    return
    raise ValueError('Failed to find Cl in attribution map')

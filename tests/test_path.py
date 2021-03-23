import unittest
from flat import command
import math

class PathTest(unittest.TestCase):
    def test_path(self):
        from bezmerizing import Path
        p = Path([
            command.moveto(10, 10),
            command.lineto(15, 25),
            command.curveto(5, 50, 15, 50, 5, 25),
            command.closepath
        ])
        transformed = p.translate(6, 6).rotate(math.pi*0.5).scale(4)
        print(transformed)
        self.assertEqual(str(transformed),
            'Path([moveto(x=64.0000, y=-64.0000), '
            'lineto(x=124.0000, y=-84.0000), '
            'curveto(x1=224.0000, y1=-44.0000, '
            'x2=224.0000, y2=-84.0000, x=124.0000, y=-44.0000), '
            'closepath])')
        
if __name__ == '__main__':
    unittest.main()

